#include "Tools/Node.h"
#include "Tools/UserInterface.h"
#include "DoubleObstacle.h"
#include "DrivingForce.h"
#include "GrainInfo.h"
#include "InterfaceEnergy.h"
#include "InterfaceMobility.h"
#include "PhaseField.h"
#include "Settings.h"
#include "VTK.h"

namespace opensim
{

using namespace std;

DoubleObstacle::DoubleObstacle(Settings& locSettings)
{
    this->Initialize(locSettings);
}

void DoubleObstacle::Initialize(Settings& locSettings)
{
    thisclassname = "DoubleObstacle";

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void DoubleObstacle::CalculatePhaseFieldIncrements(PhaseField& Phase, InterfaceEnergy& Sigma,
                                        InterfaceMobility& Mu)
{

    const double Prefactor = Pi*Pi/(Phase.Eta*Phase.Eta);
    //const double Prefactor1 = 8.0*Phase.Eta/(Pi*Pi);                      // can be used if the curvature driving force is calculated
    //double correction = 1.0 + 2.0/(Phase.iWidth*Phase.iWidth);                  // compensates for curvature effect change with the change of the interface width
    if(Phase.NucleationPresent)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
        {
            if (Phase.Fields(i,j,k).flag)
            {
                double norm_1 = 1.0/double(Phase.Fields(i,j,k).size());

                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta < Phase.Fields(i,j,k).cend(); ++beta)
                {
                    double scaleA = 1.0;
                    double scaleB = 1.0;

                    if(Phase.FieldsStatistics[alpha->index].IsNucleus())                                   //判断是否形核
                    {
                        scaleA = Phase.FieldsStatistics[alpha->index].MAXVolume/Phase.RefVolume;
                    }
                    if(Phase.FieldsStatistics[ beta->index].IsNucleus())                                    //判断是否形核
                    {
                        scaleB = Phase.FieldsStatistics[ beta->index].MAXVolume/Phase.RefVolume;
                    }

                    double scale = sqrt(scaleA * scaleB);

                    double dPhi_dt  = Sigma(i, j, k, alpha->index, beta->index)*
                                     ((alpha->laplacian + 0.5*(1.0 + scale)*Prefactor*alpha->value) -
                                      ( beta->laplacian + 0.5*(1.0 + scale)*Prefactor* beta->value));

                    if(Phase.Fields(i,j,k).size() > 2)
                    for(auto gamma = Phase.Fields(i,j,k).cbegin();
                             gamma < Phase.Fields(i,j,k).cend(); ++gamma)
                    if((gamma != alpha) && (gamma != beta))
                    {
                        dPhi_dt += (Sigma(i, j, k,  beta->index, gamma->index) -
                                    Sigma(i, j, k, alpha->index, gamma->index))*
                                   (gamma->laplacian + 0.5*(1.0 + scale)*Prefactor*gamma->value);
                    }
                    dPhi_dt *= Mu(i,j,k, alpha->index, beta->index)*norm_1*scale;//*correction;

                    Phase.FieldsDot(i,j,k).add_asym(alpha->index, beta->index, dPhi_dt);
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    else   /// No NucleationPresent
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
        {
            if (Phase.Fields(i,j,k).flag)
            {
                double norm_1 = 1.0/double(Phase.Fields(i,j,k).size());

                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta != Phase.Fields(i,j,k).cend(); ++beta)
                {
                    double dPhi_dt = Sigma(i, j, k, alpha->index, beta->index)*
                                    ((alpha->laplacian + Prefactor*alpha->value) -
                                     ( beta->laplacian + Prefactor* beta->value));

                    if(Phase.Fields(i,j,k).size() > 2)
                    for(auto gamma = Phase.Fields(i,j,k).cbegin();
                             gamma != Phase.Fields(i,j,k).cend(); ++gamma)
                    if((gamma != alpha) && (gamma != beta))
                    {
                        dPhi_dt += (Sigma(i, j, k,  beta->index, gamma->index) -
                                    Sigma(i, j, k, alpha->index, gamma->index))*
                                   (gamma->laplacian + Prefactor*gamma->value);
                    }

                    dPhi_dt *= Mu(i,j,k, alpha->index, beta->index)*norm_1;//*correction;

                    Phase.FieldsDot(i,j,k).add_asym(alpha->index, beta->index, dPhi_dt);
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

void DoubleObstacle::FixSpreading(PhaseField& Phase, BoundaryConditions& BC, double cutoff)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        for (auto alpha = Phase.Fields(i,j,k).begin();
                  alpha != Phase.Fields(i,j,k).end(); alpha++)
        {
            if(alpha->value < cutoff and Phase.FieldsStatistics[alpha->index].Stage == 0)
            {
                alpha->value = 0.0;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Phase.Finalize(BC);
}

double DoubleObstacle::PointEnergy(PhaseField& Phase, InterfaceEnergy& Sigma,
        int i, int j, int k)
{
    double energy = 0;
    double Prefactor = Phase.Eta*Phase.Eta/(Pi*Pi);

    if (Phase.Fields(i,j,k).flag)
    {
        NodeV locGradients = Phase.Gradients(i,j,k);
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta < Phase.Fields(i,j,k).cend();  ++beta)
        {
            double alphaGradX = 0.0;
            double alphaGradY = 0.0;
            double alphaGradZ = 0.0;

            double betaGradX = 0.0;
            double betaGradY = 0.0;
            double betaGradZ = 0.0;

            locGradients.get(alpha->index, alphaGradX, alphaGradY, alphaGradZ);
            locGradients.get( beta->index,  betaGradX,  betaGradY,  betaGradZ);

            energy += 4.0*Sigma(i,j,k,alpha->index, beta->index)/Phase.Eta*
                                       (alpha->value*beta->value -
                                        Prefactor*
                                        (alphaGradX*betaGradX +
                                         alphaGradY*betaGradY +
                                         alphaGradZ*betaGradZ));
        }
    }
    return energy;
}

double DoubleObstacle::Energy(PhaseField& Phase, InterfaceEnergy& Sigma)
{
    double Energy = 0.0;
    const double dx = Phase.dx;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(+:Energy))
    {
        Energy += PointEnergy(Phase, Sigma, i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    return Energy*dx*dx*dx;
}

double DoubleObstacle::AverageEnergyDensity(PhaseField& Phase,
        InterfaceEnergy& Sigma)
{
    double Energy = 0.0;
    const double Nx = Phase.Nx;
    const double Ny = Phase.Ny;
    const double Nz = Phase.Nz;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(+:Energy))
    {
        Energy += PointEnergy(Phase, Sigma, i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    return Energy/(Nx*Ny*Nz);
}

void DoubleObstacle::WriteEnergyVTK(PhaseField& Phase, InterfaceEnergy& Sigma, int tStep)
{
    stringstream buffer;
    std::vector<int> DataTypes {PDScalars};

    VTK::WriteHeader(buffer, Phase.Nx, Phase.Ny, Phase.Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        buffer << "<DataArray type = \"Float64\" Name = \"" << "InterfaceEnergy" <<
                    "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < Phase.Nz; ++k)
        for(int j = 0; j < Phase.Ny; ++j)
        for(int i = 0; i < Phase.Nx; ++i)
        {
            buffer << PointEnergy(Phase, Sigma,i,j,k) << endl;
        }
        buffer << "</DataArray>" << endl;
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Phase.Nx, Phase.Ny, Phase.Nz);
    VTK::WriteToFile(buffer, "InterfaceEnergy", tStep);
}
}// namespace opensim