

#include "Tools/CommonFunctions.h"
#include "Tools/ElasticityTensors.h"
#include "Tools/Node.h"
#include "Tools/UserInterface.h"
#include "BoundaryConditions.h"
#include "Compositions.h"
#include "FluidDynamics/D3Q27.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "Info.h"
//#include "Parallelism/BoundaryConditionsParallel.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Velocities.h"
#include <math.h>

#define RANGE 2 ///< Averaging range for appearing and disappearing obstacles

namespace opensim
{

FlowSolverLBM::FlowSolverLBM(const Settings& locSettings,
        const std::string InputFileName, size_t boundary)
{
    this->Initialize(locSettings, InputFileName, boundary);
}

void FlowSolverLBM::Initialize(const Settings& locSettings,
        const std::string InputFileName, size_t boundary)
{
    thisclassname = "FlowSolverLBM";
    //DefaultInputFileName = ProjectInputDir + "FlowInput.opi";

    Nx      = locSettings.Nx;
    Ny      = locSettings.Ny;
    Nz      = locSettings.Nz;
    dt      = locSettings.dt;
    dx      = locSettings.dx;
    dx3     = dx*dx*dx;
    dxt     = dx/dt;
    dxt2    = dx/dt/dt;
    dtx     = dt/dx;
    Nphases = locSettings.Nphases;
    N_Comp  = locSettings.Ncomp - 1;

    if(InputFileName == std::string())
    {
        this->ReadInput();
    }
    else
    {
        this->ReadInput(DefaultInputFileName);
    }

    eta = locSettings.iWidth;

    lbForceDensity.Allocate   (Nx, Ny, Nz, {N_Fluid_Comp}, boundary);
    Obstacle.Allocate         (Nx, Ny, Nz, boundary);
    ObstacleAppeared.Allocate (Nx, Ny, Nz, boundary);
    ObstacleVanished.Allocate (Nx, Ny, Nz, boundary);
    lbDensity.Allocate        (Nx, Ny, Nz, {N_Fluid_Comp}, boundary);
    lbMomentum.Allocate       (Nx, Ny, Nz, {N_Fluid_Comp}, boundary);
    lbMomentumDelta.Allocate  (Nx, Ny, Nz, {N_Fluid_Comp}, boundary);
    lbPopulations.Allocate    (Nx, Ny, Nz, {N_Fluid_Comp}, boundary);
    lbPopulationsTMP.Allocate (Nx, Ny, Nz, {N_Fluid_Comp}, boundary);

    FluidMass.resize(N_Fluid_Comp);
    //V_star.resize   (N_Fluid_Comp); //Note Unused
    tau.resize      (N_Fluid_Comp);

    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        FluidMass    [n] = 0;
        //V_star       [n] = 36.0*nu[n]*h_star/(eta*eta);
        tau          [n] = nu[n]*dt*3.0/(dx*dx) + 0.5;
        if(tau[n] < 0.53)
        {
            Info::WriteWarning("tau is too small!", thisclassname, "Initialize");
        }
    }
    //std::cout << "G[0][0] = " << G[0][0]  << std::endl; //NOTE is done in ReadInput
    //std::cout << "rho_0 = "   << rho_0[0] << std::endl; //NOTE is done in ReadInput

    Info::WriteStandard("tau", tau[0]);

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void FlowSolverLBM::Read(const int tStep)
{
    std::string FileName = UserInterface::MakeFileName(RawDataDir,"LBM_", tStep, ".dat");

    std::fstream inp(FileName.c_str(), std::ios::in | std::ios::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be opened", thisclassname, "Read");
        exit(1);
    };

    STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,lbPopulations.Bcells())
    {
        for(int n = 0; n < N_Fluid_Comp; ++n)
        for(int ii = -1; ii <= 1; ++ii)
        for(int jj = -1; jj <= 1; ++jj)
        for(int kk = -1; kk <= 1; ++kk)
        {
            inp.read(reinterpret_cast<char*>(&lbPopulations(i,j,k)({n})(ii,jj,kk)), sizeof(double));
        }
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,Obstacle,Obstacle.Bcells())
    {
        inp.read(reinterpret_cast<char*>(&Obstacle(i,j,k)), sizeof(bool));
    }
    STORAGE_LOOP_END

    CalcateDensityAndMomentum();

    Info::WriteStandard(thisclassname, "Binary input loaded");
}

void FlowSolverLBM::ReadInput(const std::string InputFileName)
{
    Info::WriteBlankLine();
    Info::WriteLineInsert("Flow solver LBM");
    Info::WriteStandard("Source", InputFileName);

    std::fstream inp(InputFileName.c_str(), std::ios::in);
    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", "ReadInput");
        exit(1);
    };

    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    for (int i = 0; i < 3; i++)
    {
        std::stringstream GAstr;
        GAstr << "G0[" << i << "]";
        GA[i] = UserInterface::ReadParameterD(inp, moduleLocation, GAstr.str());
    }
    N_Fluid_Comp = UserInterface::ReadParameterD(inp, moduleLocation, std::string("N_FLUID_COMP"));
    h_star = UserInterface::ReadParameterD(inp, moduleLocation, std::string("H_STAR"));
    G.resize(N_Fluid_Comp);
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        G[n].resize(N_Fluid_Comp);
        for (int m = 0; m < N_Fluid_Comp; ++m)
        {
            std::stringstream GBstr;
            GBstr << "GB[" << n << "][" << m << "]";
            G[n][m] = UserInterface::ReadParameterD(inp, moduleLocation, GBstr.str());
        }
    }
    nu.resize(N_Fluid_Comp);
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        std::stringstream NUstr;
        NUstr << "NU[" << n << "]";
        nu[n] = UserInterface::ReadParameterD(inp, moduleLocation, NUstr.str());
    }
    rho_0.resize(N_Fluid_Comp);
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        std::stringstream RHO0str;
        RHO0str << "RHO0[" << n << "]";
        rho_0[n] = UserInterface::ReadParameterD(inp, moduleLocation, RHO0str.str(), false, 1.0);
    }

    if (N_Comp > 0)
    {
        drhodc.resize(N_Comp);
        for (int n = 0; n < N_Comp; ++n)
        {
            std::stringstream DRHODCstr;
            DRHODCstr << "DRHODC[" << n << "]";
            drhodc[n] = UserInterface::ReadParameterD(inp, moduleLocation, DRHODCstr.str());
        }
    }

    //RhoSol = UserInterface::ReadParameterD(inp, moduleLocation, std::string("RhoSol"),false, 0.0);  //UNUSED
    //RhoLiq = UserInterface::ReadParameterD(inp, moduleLocation, std::string("RhoLiq"),false, 0.0);  //UNUSED
    dRho = UserInterface::ReadParameterD(inp, moduleLocation, std::string("dRho"));

    Do_Benzi       = UserInterface::ReadParameterB(inp, moduleLocation, std::string("BENZI"));
    Do_Gravitation = UserInterface::ReadParameterB(inp, moduleLocation, std::string("GRAVITY"));
    Do_Buoyancy    = UserInterface::ReadParameterB(inp, moduleLocation, std::string("BUOYANCY"));
    Do_Drag        = UserInterface::ReadParameterB(inp, moduleLocation, std::string("DRAG"));
    Do_BounceBack  = UserInterface::ReadParameterB(inp, moduleLocation, std::string("BOUNCEBACK"));
    Do_SolidSolid  = UserInterface::ReadParameterB(inp, moduleLocation, std::string("SOLIDSOLID"));
    inp.close();
    lbGA = GA*dt*dtx;

    Info::WriteLine();
}

void FlowSolverLBM::SetUniformVelocity(const BoundaryConditions& BC,
        const dVector3 U0)
{
    dVector3 lbU0 = U0*dtx;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,lbDensity.Bcells(),)
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        lbDensity(i,j,k)({n}) = 1.0;
        lbMomentum(i,j,k)({n}) = lbU0;
        lbPopulations(i,j,k)({n}) = D3Q27::EqDistribution(1.0, lbMomentum(i,j,k)({n}));
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    SetBoundaryConditions(BC);
}

void FlowSolverLBM::Write(const int tStep) const
{
    std::string FileName = UserInterface::MakeFileName(RawDataDir,"LBM_",tStep,".dat");

    std::ofstream out(FileName.c_str(), std::ios::out | std::ios::binary);

    if (!out)
    {
        std::cout << "File \"" << FileName << "\" could not be created! Terminating!!!" << std::endl;
        exit(1);
    };

    STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,lbPopulations.Bcells())
    {
        for(int n = 0; n < N_Fluid_Comp; ++n)
        for(int ii = -1; ii <= 1; ++ii)
        for(int jj = -1; jj <= 1; ++jj)
        for(int kk = -1; kk <= 1; ++kk)
        {
            double value = lbPopulations(i,j,k)({n})(ii,jj,kk);
            out.write(reinterpret_cast<char*>(&value), sizeof(double));
        }
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,Obstacle,Obstacle.Bcells())
    {
        bool value = Obstacle(i,j,k);
        out.write(reinterpret_cast<char*>(&value), sizeof(bool));
    }

    STORAGE_LOOP_END

    out.close();
}

void FlowSolverLBM::WriteVTK(const int tStep) const
{
    std::stringstream outbufer;
    outbufer.precision(5);
    outbufer << "# vtk DataFile Version 3.0\n";
    outbufer << "PhaseField\n";
    outbufer << "ASCII\n";
    outbufer << "DATASET STRUCTURED_GRID\n";
    outbufer << "DIMENSIONS " << Nx << " " << Ny  << " " << Nz  << "\n";
    outbufer << "POINTS " << Nx*Ny*Nz << " int\n";
    for (int k = 0; k < Nz; ++k)
    for (int j = 0; j < Ny; ++j)
    for (int i = 0; i < Nx; ++i)
    {
        outbufer << i << " " << j  << " " << k  << "\n";
    }
    outbufer << "\n";
    outbufer << "POINT_DATA " << Nx*Ny*Nz << "\n";
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        outbufer << "SCALARS Density_" << n << " double" << std::endl;
        outbufer << "LOOKUP_TABLE default" << "\n";
        for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx; ++i)
        {
            outbufer <<  lbDensity(i, j, k)({n}) << "\n";
        }
        outbufer << "\n";
    }
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        outbufer << "SCALARS VX_" << n << " double" << std::endl;
        outbufer << "LOOKUP_TABLE default" << "\n";
        for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx; ++i)
        {
            if (lbDensity(i,j,k)({n}) > DBL_EPSILON )
            {
                outbufer << (lbMomentum(i, j, k)({n})[0]
                        + 0.5*lbMomentumDelta(i, j, k)({n})[0])/lbDensity(i, j, k)({n}) << "\n";
            }
            else
            {
                outbufer << 0.0 << "\n";
            }
        }
        outbufer << "\n";
    }
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        outbufer << "SCALARS VY_" << n << " double" << std::endl;
        outbufer << "LOOKUP_TABLE default" << "\n";
        for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx; ++i)
        {
            if (lbDensity(i,j,k)({n}) > DBL_EPSILON )
            {
                outbufer << (lbMomentum(i, j, k)({n})[1]
                        + 0.5*lbMomentumDelta(i, j, k)({n})[1])/lbDensity(i, j, k)({n}) << "\n";
            }
            else
            {
                outbufer << 0.0 << "\n";
            }
        }
        outbufer << "\n";
    }
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        outbufer << "SCALARS VZ_" << n << " double" << std::endl;
        outbufer << "LOOKUP_TABLE default" << "\n";
        for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx; ++i)
        {
            if (lbDensity(i,j,k)({n}) > DBL_EPSILON )
            {
                outbufer << (lbMomentum(i, j, k)({n})[2]
                        + 0.5*lbMomentumDelta(i, j, k)({n})[2])/lbDensity(i, j, k)({n}) << "\n";
            }
            else
            {
                outbufer << 0.0 << "\n";
            }
        }
        outbufer << "\n";
    }
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        outbufer << "SCALARS Pressure" << n << " double" << std::endl;
        outbufer << "LOOKUP_TABLE default" << "\n";
        for (int k = 0; k < Nz; ++k)
        for (int j = 0; j < Ny; ++j)
        for (int i = 0; i < Nx; ++i)
        {
            if (lbDensity(i, j, k)({n}) > DBL_EPSILON)
            {
                outbufer << lbDensity(i, j, k)({n}) / 3. + 0.5 * G[n][0] *
                    psi(lbDensity(i,j,k)({n}),rho_0[n]) *
                    psi(lbDensity(i,j,k)({n}),rho_0[n]) / 3. << "\n";
            }
            else
            {
                outbufer << 0.0 << "\n";
            }
        }
        outbufer << "\n";
    }

    outbufer << "\n";
    outbufer << "SCALARS Obstacle double 1" << std::endl;
    outbufer << "LOOKUP_TABLE default" << "\n";
    for (int k = 0; k < Nz; ++k)
    for (int j = 0; j < Ny; ++j)
    for (int i = 0; i < Nx; ++i)
    {
        outbufer << Obstacle(i,j,k) << "\n";
    }
    outbufer << "\n";
    std::ostringstream FileName;
    FileName.fill('0');
    FileName << VTKDir << "LBDensity" << std::setw(8) << tStep << ".vtk";

    std::ofstream vtk_file(FileName.str());
    vtk_file << outbufer.rdbuf();
    vtk_file.close();
}

void FlowSolverLBM::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetXVector(lbPopulations);
    BC.SetYVector(lbPopulations);
    BC.SetZVector(lbPopulations);
    BC.SetX(lbDensity);
    BC.SetY(lbDensity);
    BC.SetZ(lbDensity);
    BC.SetXVector(lbMomentum);
    BC.SetYVector(lbMomentum);
    BC.SetZVector(lbMomentum);
}

void FlowSolverLBM::CalculateForceBenzi(const int i,const int j,const int k,
        PhaseField& Phase)
{
    for (int ii = -1; ii < 2; ++ii)
    for (int jj = -1; jj < 2; ++jj)
    for (int kk = -1; kk < 2; ++kk)
    for (int n = 0; n < N_Fluid_Comp; ++n)
    for (int m = 0; m < N_Fluid_Comp; ++m)
    {
        const double f = - G[n][m] * psi(lbDensity(i,j,k)({n}), rho_0[n]) *
            (D3Q27::weights[ii+1][jj+1][kk+1] *
             (psi(lbDensity(i+ii,j+jj,k+kk)({m}), rho_0[m])));

        const dVector3 lbBenziForceDensity = {f*ii,f*jj,f*kk};

        lbForceDensity(i,j,k)({n}) += lbBenziForceDensity/dt;

        if (!Obstacle(i,j,k) && Obstacle(i+ii,j+jj,k+kk))
        for(auto alpha = Phase.Fields(i+ii,j+jj,k+kk).cbegin();
             alpha != Phase.Fields(i+ii,j+jj,k+kk).cend(); ++alpha)
        if(Phase.FieldsStatistics[alpha->index].State == Solid &&
           Phase.FieldsStatistics[alpha->index].IsMobile)
        {
            const dVector3 pos = {double(i+ii), double(j+jj), double(k+kk)};
            dVector3 distanceCM;
            CommonFunctions::CalculateDistancePeriodic(pos,
                    Phase.FieldsStatistics[alpha->index].Rcm, distanceCM,
                    Nx, Ny, Nz);
            const dVector3 R = distanceCM * dx;
            const dVector3 BenziForceDensity = lbBenziForceDensity*dRho*dxt2*(-1);

            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                Phase.FieldsStatistics[alpha->index].Force
                    += BenziForceDensity * dx3;
                Phase.FieldsStatistics[alpha->index].Torque
                    += R.cross(BenziForceDensity) * dx3;
            }
        }
    }
}

void FlowSolverLBM::CalculateForceGravitation(const int i, const int j,
        const int k)
{
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        lbForceDensity(i,j,k)({n})[0] += lbGA[0]*lbDensity(i,j,k)({n})/dt;
        lbForceDensity(i,j,k)({n})[1] += lbGA[1]*lbDensity(i,j,k)({n})/dt;
        lbForceDensity(i,j,k)({n})[2] += lbGA[2]*lbDensity(i,j,k)({n})/dt;
    }
}

void FlowSolverLBM::CalculateForceBuoyancy(const int i, const int j,
        const int k, PhaseField& Phase, const Composition& Cx)
{
    for(auto it = Phase.Fields(i,j,k).cbegin();
             it != Phase.Fields(i,j,k).cend(); ++it)
    {
        int pIndex = Phase.FieldsStatistics[it->index].Phase;
        if(Phase.FieldsStatistics[it->index].State == Liquid)
        for (int n = 0; n < N_Fluid_Comp; ++n)
        for (int Comp = 0; Comp < N_Comp; ++Comp)
        {
            lbForceDensity(i,j,k)({n}) += lbGA * lbDensity(i,j,k)({n}) *
                drhodc[Comp]*(Cx.Phase(i,j,k)({pIndex, Comp}) -
                        Cx.Initial({pIndex, Comp}))/dt;
        }
    }
}

void FlowSolverLBM::CalculateForceDrag(const int i,const int j,const int k,
        PhaseField& Phase, const Velocities& Vel)
{
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        //Solved for drag Force:
        if(lbDensity(i,j,k)({n}) != 0.0)
        {
            double SolidFraction = 0.0;
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            if(Phase.FieldsStatistics[alpha->index].State == Solid)
            {
                SolidFraction += alpha->value;
            }
            double loclbDensity_1 = 1.0/lbDensity(i,j,k)({n});
            double mu = nu[n]*lbDensity(i,j,k)({n});  //Dynamic viscosity
            for(auto it = Phase.Fields(i,j,k).cbegin();
                     it != Phase.Fields(i,j,k).cend(); ++it)
            {
                if(Phase.FieldsStatistics[it->index].State == Solid)
                {
                    const double factor =
                        (it->value * it->value * (1.0 - SolidFraction))/
                        (Phase.iWidth * Phase.iWidth);
                    const dVector3 lbDragForce =
                        (lbMomentum(i,j,k)({n}) * loclbDensity_1 -
                         Vel.Phase(i,j,k)({Phase.FieldsStatistics[it->index].Phase}) *
                         dtx ) * mu * factor * h_star;
                    lbForceDensity(i,j,k)({n}) -= lbDragForce;

                    dVector3 locR;
                    locR[0] = (i - Phase.FieldsStatistics[it->index].Rcm[0])*dx;
                    locR[1] = (j - Phase.FieldsStatistics[it->index].Rcm[1])*dx;
                    locR[2] = (k - Phase.FieldsStatistics[it->index].Rcm[2])*dx;
                    #ifdef _OPENMP
                    #pragma omp critical
                    #endif
                    {
                        Phase.FieldsStatistics[it->index].Force  += lbDragForce * dRho * dxt2 * dx3;
                        Phase.FieldsStatistics[it->index].Torque += locR.cross(lbDragForce) * dRho * dxt2 * dx3;
                    }
                }
            }
        }
    }
}

double FlowSolverLBM::BounceBack(const int i, const int j, const int k,
        const int ii, const int jj, const int kk, const int n,
        PhaseField& Phase, const BoundaryConditions& BC, double& lbDensityChange)
{
    double lbG = 0.;
    if (Do_BounceBack)
    for(auto alpha  = Phase.Fields(i-ii,j-jj,k-kk).cbegin();
             alpha != Phase.Fields(i-ii,j-jj,k-kk).cend(); ++alpha)
    if(Phase.FieldsStatistics[alpha->index].State == Solid &&
       Phase.FieldsStatistics[alpha->index].IsMobile)
    {
        dMatrix3x3 W;
        W(0,0) = 0.0;
        W(1,1) = 0.0;
        W(2,2) = 0.0;
        W(0,1) = -Phase.FieldsStatistics[alpha->index].aVel[2];
        W(0,2) =  Phase.FieldsStatistics[alpha->index].aVel[1];
        W(1,2) = -Phase.FieldsStatistics[alpha->index].aVel[0];
        W(1,0) =  Phase.FieldsStatistics[alpha->index].aVel[2];
        W(2,0) = -Phase.FieldsStatistics[alpha->index].aVel[1];
        W(2,1) =  Phase.FieldsStatistics[alpha->index].aVel[0];

        const dVector3 pos = {i-0.5*ii, j-0.5*jj, k-0.5*kk};
        dVector3 distanceCM;
        CommonFunctions::CalculateDistancePeriodic(pos,
                Phase.FieldsStatistics[alpha->index].Rcm,distanceCM, Nx, Ny, Nz);
        const dVector3 R   = distanceCM * dx;
        const dVector3 Vel = Phase.FieldsStatistics[alpha->index].Vcm + W * R;

        const double tmp = D3Q27::weights[ii+1][jj+1][kk+1] *
            lbDensity(i,j,k)({n}) * (Vel[0]*ii + Vel[1]*jj + Vel[2]*kk) * dtx;

        const double lbBBDensity =
            2.0 * (lbPopulations(i,j,k)({n})(-ii,-jj,-kk) + 3.0 * tmp);

        lbG += alpha->value * 6.0 * tmp;
        //lbG += 6.0 * tmp;  // Unstable

        // NOTE: Bounce Back Momentum Density is now in physical units
        dVector3 BounceBackForceDensity;
        BounceBackForceDensity[0] = - ii * lbBBDensity * dRho * dxt2;
        BounceBackForceDensity[1] = - jj * lbBBDensity * dRho * dxt2;
        BounceBackForceDensity[2] = - kk * lbBBDensity * dRho * dxt2;

        #ifdef _OPENMP
        #pragma omp critical
        #endif
        {
            Phase.FieldsStatistics[alpha->index].Force  +=
                BounceBackForceDensity * dx3;

            Phase.FieldsStatistics[alpha->index].Torque +=
                R.cross(BounceBackForceDensity * dx3);
        }
    }
    lbDensityChange += lbG;

    // Return bounced back population from resting solid boundary and the
    // correction according to the movement of the solid boundary.
    return lbPopulations(i,j,k)({n})(-ii,-jj,-kk) + lbG;
}

void FlowSolverLBM::Propagation(PhaseField& Phase, const BoundaryConditions& BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,0,)
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        double lbDensityChange = 0;
        for (int ii = -1; ii < 2; ii++)
        for (int jj = -1; jj < 2; jj++)
        for (int kk = -1; kk < 2; kk++)
        {
            if (Obstacle(i, j, k))
            {
                lbPopulationsTMP(i,j,k)({n})(ii,jj,kk) =
                    lbPopulations(i,j,k)({n})(ii,jj,kk);
            }
            else if (Obstacle(i-ii, j-jj, k-kk))
            {
                lbPopulationsTMP(i,j,k)({n})(ii,jj,kk) =
                    BounceBack(i, j, k, ii, jj, kk, n, Phase, BC, lbDensityChange);
            }
            else
            {
                lbPopulationsTMP(i,j,k)({n})(ii,jj,kk) =
                    lbPopulations(i-ii, j-jj, k-kk)({n})(ii,jj,kk);
            }
        }
        lbPopulationsTMP(i,j,k)({n})(0,0,0) -= lbDensityChange; //NoSlip
        //lbPopulationsTMP(i,j,k)({n}) +=
        //    D3Q27::EqDistribution(lbDensity(i,j,k)({n})) -
        //    D3Q27::EqDistribution(lbDensity(i,j,k)({n}) + lbDensityChange);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,0,)
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        lbPopulations(i,j,k)({n}) = lbPopulationsTMP(i,j,k)({n});
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::ApplyForces(PhaseField& Phase, const Composition& Cx,
        const Velocities& Vel, const BoundaryConditions& BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,0,)
    if (!Obstacle(i,j,k))
    {
        for (int n = 0; n < N_Fluid_Comp; ++n)
        {
            lbForceDensity(i,j,k)({n}).set_to_zero();
        }
        if (Do_Benzi) CalculateForceBenzi(i,j,k, Phase);
        if (Do_Gravitation) CalculateForceGravitation(i,j,k);
        if (Do_Buoyancy) CalculateForceBuoyancy(i,j,k, Phase, Cx);
        if (Do_Drag and Phase.Interface(i,j,k))
        {
            CalculateForceDrag(i,j,k, Phase, Vel);
        }
        for (int n = 0; n < N_Fluid_Comp; ++n)
        {
            lbMomentumDelta(i,j,k)({n}) = lbForceDensity(i,j,k)({n})*dt;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::ApplyForces(PhaseField& Phase, const Velocities& Vel,
        const BoundaryConditions& BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,0,)
    if (!Obstacle(i,j,k))
    {
        for (int n = 0; n < N_Fluid_Comp; ++n)
        {
            lbForceDensity(i,j,k)({n}).set_to_zero();
        }
        if (Do_Benzi) CalculateForceBenzi(i,j,k, Phase);
        if (Do_Gravitation) CalculateForceGravitation(i,j,k);
        if (Do_Drag and Phase.Interface(i,j,k))
        {
            CalculateForceDrag(i,j,k, Phase, Vel);
        }
        for (int n = 0; n < N_Fluid_Comp; ++n)
        {
            lbMomentumDelta(i,j,k)({n}) = lbForceDensity(i,j,k)({n})*dt;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::Collision()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,0,)
    if (!Obstacle(i,j,k))
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
       const double q = 1.0 - 1.0/tau[n];
       lbPopulations(i,j,k)({n}) = (lbPopulations(i,j,k)({n}) -
                                    D3Q27::EqDistribution(lbDensity(i,j,k)({n}), lbMomentum(i,j,k)({n})))*q +
                                    D3Q27::EqDistribution(lbDensity(i,j,k)({n}), (lbMomentum(i,j,k)({n}) + lbMomentumDelta(i,j,k)({n})));
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

size_t FlowSolverLBM::CountObstacleNodes(void) const
{
    size_t ObstacleNodes = 0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,0, reduction(+:ObstacleNodes))
    {
        if (Obstacle(i,j,k))
        {
            ObstacleNodes++;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    return ObstacleNodes;
}

size_t FlowSolverLBM::CountFluidNodes(void) const
{
    size_t FluidNodes = 0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,0, reduction(+:FluidNodes))
    {
        if (not Obstacle(i,j,k))
        {
            FluidNodes++;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    return FluidNodes;
}

void FlowSolverLBM::EnforceMassConservation(void)
{
    // Calculate Fluid mass
    std::vector<double> DeltaDensity = CalculateFluidMass();
    FluidNodes = CountFluidNodes();

    for (int n = 0; n < N_Fluid_Comp; ++n)
    if (FluidMass[n] != 0)
    {
        DeltaDensity[n] -= FluidMass[n];
        DeltaDensity[n] /= FluidNodes;

        // Fix density
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,0,)
        {
            if (not Obstacle(i,j,k))
            {
                lbPopulations(i,j,k)({n}) = lbPopulations(i,j,k)({n}) -
                    D3Q27::EqDistribution(lbDensity(i,j,k)({n}) + DeltaDensity[n]) +
                    D3Q27::EqDistribution(lbDensity(i,j,k)({n}));
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    FluidMass = CalculateFluidMass();
}

void FlowSolverLBM::EnforceSolidMomentum(PhaseField& Phase, const dVector3 value)
{
    // Fix Solid momtum
    const double dx   = Phase.dx;
    const size_t size = Phase.FieldsStatistics.size();

    dVector3 SolidMomentum = {0.0,0.0,0.0};
    std::vector<double> SolidMass(size,0.0);
    size_t NumSolids = 0;
    for (auto grain : Phase.FieldsStatistics.GrainStorage)
    if (grain.State == Solid and grain.IsMobile)
    {
        SolidMass[NumSolids] = grain.Volume * grain.Density*dx*dx*dx;
        SolidMomentum       += grain.Vcm * SolidMass[NumSolids];
        NumSolids++;
    }

    const dVector3 SolidMomentumFix = (SolidMomentum - value)/NumSolids;
    NumSolids = 0;
    for (auto& grain : Phase.FieldsStatistics.GrainStorage)
    if (grain.State == Solid and grain.IsMobile)
    {
        grain.Vcm -= SolidMomentumFix/SolidMass[NumSolids];
        NumSolids++;
    }
}

void FlowSolverLBM::DetectObstaclesSimple(const PhaseField& Phase)
{
    ObstaclesChanged = false;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,lbDensity.Bcells(),)
    {
        // Detect the number of present solids
        for(auto it = Phase.Fields(i,j,k).cbegin();
                it != Phase.Fields(i,j,k).cend(); ++it)
        if(Phase.FieldsStatistics[it->index].State == Solid)
        {
            if(it->value >= 0.95)
            {
                Obstacle(i,j,k) = true;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::DetectObstacles(const PhaseField& Phase)
{
    // Check if PhaseField has enough boundary cells
    if (Phase.Fields.Bcells() < 2*RANGE)
    {
        const std::string message = "Allocate at least " +
            std::to_string(2*RANGE) +
            " boundary cells for the phase-field!";

        Info::WriteExit(message, thisclassname, "DetectObstacles");
    }

    ObstaclesChanged = false;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,lbDensity.Bcells(),)
    {
        ObstacleVanished(i,j,k) = false;
        ObstacleAppeared(i,j,k) = false;

        // Detect the number of present solids
        //double SolidFraction = 0.0;
        size_t num_solids = 0;
        for(auto it = Phase.Fields(i,j,k).cbegin();
                it != Phase.Fields(i,j,k).cend(); ++it)
        if(Phase.FieldsStatistics[it->index].State == Solid)
        {
            //SolidFraction += it-value;
            if(it->value >= 0.95)
            {
                num_solids++;
            }
        }

        // Set the new obstacles and record change
        //if(SolidFraction >= 0.95)
        if(num_solids >= 1)
        {
            if(not Obstacle(i,j,k))
            {
                ObstacleAppeared(i,j,k) = true;
                ObstaclesChanged = true;
            }
            Obstacle(i,j,k) = true;
        }
        else
        {
            if(Obstacle(i,j,k))
            {
                ObstacleVanished(i,j,k) = true;
                ObstaclesChanged = true;
            }
            Obstacle(i,j,k) = false;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

dVector3 FlowSolverLBM::CalculateFluidMomentum(const PhaseField& Phase) const
{
    dVector3 Momentum = {0.,0.,0.};

    // Calculate fluid momentum
    #pragma omp parallel
    {
        dVector3 locMomentum = {0.,0.,0.};
        #pragma omp for collapse(OMP_COLLAPSE_LOOPS)\
        schedule(dynamic,OMP_DYNAMIC_CHUNKSIZE)
        for (int i = 0; i < Nx; i++)
        for (int j = 0; j < Ny; j++)
        for (int k = 0; k < Nz; k++)
        if (not Obstacle(i,j,k))
        for (int n = 0; n < N_Fluid_Comp; ++n)
        for (int ii = -1; ii < 2; ii++)
        for (int jj = -1; jj < 2; jj++)
        for (int kk = -1; kk < 2; kk++)
        {
            double rr = lbPopulations(i,j,k)({n})(ii,jj,kk);
            locMomentum[0] += rr*ii;
            locMomentum[1] += rr*jj;
            locMomentum[2] += rr*kk;
        }
        #pragma omp critical
        {
            Momentum += locMomentum;
        }
    }

    return Momentum * dx /dt *dRho*dx*dx*dx;
}

void FlowSolverLBM::SetFluidNodesNearObstacle()
{
    if(ObstaclesChanged)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,lbPopulations.Bcells()-RANGE,)
        for (int n = 0; n < N_Fluid_Comp; ++n)
        {
            lbPopulationsTMP(i,j,k)({n}).set_to_zero();
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,lbPopulations.Bcells()-RANGE,)
        if (ObstacleVanished(i,j,k))
        for (int n = 0; n < N_Fluid_Comp; ++n)
        {
            // If an obstacle vanishes the present cell has to be filled with the
            // surrounding fluid form the neighbour cells

            // Calculate local average population
            double   avDensity  = 0;
            dVector3 avMomentum = {0.0,0.0,0.0};
            int FluidNeighbours = 0;
            for (long int ii = -RANGE; ii <= RANGE; ++ii)
            for (long int jj = -RANGE; jj <= RANGE; ++jj)
            for (long int kk = -RANGE; kk <= RANGE; ++kk)
            if (not Obstacle(i+ii,j+jj,k+kk)
                    and not ObstacleAppeared(i+ii,j+jj,k+kk)
                    and not ObstacleVanished(i+ii,j+jj,k+kk)
#if (RANGE > 1)
                    and ii*ii + jj*jj + kk*kk <= RANGE*RANGE
#endif
                    and ii + jj + kk != 0)
            {
                avDensity  += lbDensity(i+ii,j+jj,k+kk)({n});
                avMomentum += lbMomentum(i+ii,j+jj,k+kk)({n});
                FluidNeighbours++;
            }
            if (FluidNeighbours)
            {
                // Calculate average population (+1 in the denominator
                // ensures that the new fluid node is equal to surrounding
                // density even, if the there is only one fluid neighbour!)
                avDensity  /=(FluidNeighbours+1);
                avMomentum /=(FluidNeighbours+1);

                const D3Q27 avPopulation =
                    D3Q27::EqDistribution(avDensity, avMomentum);

                // Set (i,j,k) to local average population
                lbPopulations(i,j,k)({n}) = avPopulation;

                // Subtract added fluid from neighbouring cells
                D3Q27 subsPopulation = avPopulation / FluidNeighbours;
                for (long int ii = -RANGE; ii <= RANGE; ++ii)
                for (long int jj = -RANGE; jj <= RANGE; ++jj)
                for (long int kk = -RANGE; kk <= RANGE; ++kk)
                if (not Obstacle(i+ii,j+jj,k+kk)
                        and not ObstacleAppeared(i+ii,j+jj,k+kk)
                        and not ObstacleVanished(i+ii,j+jj,k+kk)
#if (RANGE > 1)
                        and ii*ii + jj*jj + kk*kk <= RANGE*RANGE
#endif
                        and ii + jj + kk != 0)
                {
                   #pragma omp critical
                   {
                       lbPopulationsTMP(i+ii,j+jj,k+kk)({n}) -= subsPopulation;
                   }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,lbPopulations.Bcells()-RANGE,)
        if (ObstacleAppeared(i,j,k))
        for (int n = 0; n < N_Fluid_Comp; ++n)
        {
        // If an obstacle appears the present fluid as to be moved to the
        // neighbouring cell

            // Count number of fluid neighbours
            int FluidNeighbours = 0;
            for (long int ii = -RANGE; ii <= RANGE; ++ii)
            for (long int jj = -RANGE; jj <= RANGE; ++jj)
            for (long int kk = -RANGE; kk <= RANGE; ++kk)
            if (not Obstacle(i+ii,j+jj,k+kk)
                    and not ObstacleAppeared(i+ii,j+jj,k+kk)
#if (RANGE > 1)
                    and ii*ii + jj*jj + kk*kk <= RANGE*RANGE
#endif
                    and ii + jj + kk != 0)
            {
                FluidNeighbours++;
            }

            if (FluidNeighbours)
            {
                // Calculate average density
                double   addDensity  = lbDensity(i,j,k)({n})/FluidNeighbours;
                dVector3 addMomentum = lbMomentum(i,j,k)({n})/FluidNeighbours;

                // Calculate fluid population which has to be added to
                // neighbour cells
                const D3Q27 addPopulation =
                    D3Q27::EqDistribution(addDensity, addMomentum);

                // Add fluid from this cell to the neighbouring cells
                for (long int ii = -RANGE; ii <= RANGE; ++ii)
                for (long int jj = -RANGE; jj <= RANGE; ++jj)
                for (long int kk = -RANGE; kk <= RANGE; ++kk)
                if (not Obstacle(i+ii,j+jj,k+kk)
                        and not ObstacleAppeared(i+ii,j+jj,k+kk)
#if (RANGE > 1)
                       and ii*ii + jj*jj + kk*kk <= RANGE*RANGE
#endif
                       and ii + jj + kk != 0)
               {
                   #pragma omp critical
                   {
                       lbPopulationsTMP(i+ii,j+jj,k+kk)({n}) += addPopulation;
                   }
               }
           }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        // Add change in fluid nodes
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,0,)
        {
            if (not Obstacle(i,j,k))
            for (int n = 0; n < N_Fluid_Comp; ++n)
            {
                // Apply fluid mass fix
                lbPopulations(i,j,k)({n}) += lbPopulationsTMP(i,j,k)({n});
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

void FlowSolverLBM::FixPopulations(void)
{
    CalcateDensityAndMomentum();

    // Add change in fluid nodes
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbPopulations,0,)
    {
        if (not Obstacle(i,j,k))
        for (int n = 0; n < N_Fluid_Comp; ++n)
        {
            bool FixPopulation = false;
            // Check for negative populations and fix them
            for (long int ii = -1; ii <= 0; ++ii)
            for (long int jj = -1; jj <= 0; ++jj)
            for (long int kk = -1; kk <= 0; ++kk)
            if (ii + jj + kk != 0)
            if (lbPopulations(i,j,k)({n})(ii,jj,kk) <= 0.0)
            {
                FixPopulation = true;
            }

            if (lbPopulations(i,j,k)({n})(0,0,0) <= 0.0)
            {
                FixPopulation = true;
            }

            if (FixPopulation)
            {
                lbPopulations(i,j,k)({n}) =
                    D3Q27::EqDistribution(lbDensity(i,j,k)({n}), lbMomentum(i,j,k)({n}));
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::SetObstacleNodes(const PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,0,)
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        if (Obstacle(i,j,k))
        {
            if (Do_Benzi)
            {
                lbDensity(i,j,k)({n}) = 0.0;
                lbMomentum(i,j,k)({n}).set_to_zero();
                double SolidFraction = 0.0;
                for (auto alpha = Phase.Fields(i,j,k).cbegin();
                          alpha != Phase.Fields(i,j,k).cend(); ++alpha)
                {
                    if (Phase.FieldsStatistics[alpha->index].State == Solid)
                    {
                        lbDensity(i,j,k)({n}) += alpha->value
                            * Phase.FieldsStatistics[alpha->index].WettingParameter;
                        // NOTE: Commented code because it is not used
                        //if (Phase.FieldsStatistics[alpha->index].IsMobile)
                        //{
                        //    const int idx = Phase.FieldsStatistics[alpha->index].Phase;
                        //    lbMomentum(i,j,k)({n}) += Vel.Phase(i,j,k)({idx}) * dtx
                        //        * alpha->value * lbDensity(i,j,k)({n});
                        //}
                        SolidFraction += alpha->value;
                    }
                }
                if (SolidFraction)
                {
                    lbDensity(i,j,k)({n})  /= SolidFraction;
                    //lbMomentum(i,j,k)({n}) /= SolidFraction;
                }
            }
            else
            {
                lbDensity(i,j,k)({n}) = 1.0;
                lbMomentum(i,j,k)({n}).set_to_zero();
            }
            lbPopulations(i,j,k)({n}) =
                D3Q27::EqDistribution(lbDensity(i,j,k)({n}), lbMomentum(i,j,k)({n}));
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::Solve(PhaseField& Phase, const Composition& Cx,
        Velocities& Vel, const BoundaryConditions& BC)
{
    for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        Phase.FieldsStatistics[idx].Force.set_to_zero();
        Phase.FieldsStatistics[idx].Torque.set_to_zero();
    }

    FixPopulations();

    DetectObstaclesSimple(Phase);
    SetObstacleNodes(Phase);
    SetBoundaryConditions(BC);

    Propagation(Phase, BC);
    CalcateDensityAndMomentum();
    SetBoundaryConditions(BC);
    ApplyForces(Phase, Cx, Vel, BC);
    Collision();

    SetBoundaryConditions(BC);

    CalculateFluidVelocities(Phase, Vel);
}

void FlowSolverLBM::Solve(PhaseField& Phase, Velocities& Vel,
        const BoundaryConditions& BC)
{
    for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        Phase.FieldsStatistics[idx].Force.set_to_zero();
        Phase.FieldsStatistics[idx].Torque.set_to_zero();
    }

    FixPopulations();

    DetectObstaclesSimple(Phase);
    SetObstacleNodes(Phase);
    SetBoundaryConditions(BC);

    Propagation(Phase,BC);
    CalcateDensityAndMomentum();
    SetBoundaryConditions(BC);
    ApplyForces(Phase, Vel, BC);
    Collision();

    SetBoundaryConditions(BC);

    CalculateFluidVelocities(Phase, Vel);
}

std::vector<double> FlowSolverLBM::CalculateFluidMass(void) const
{
    std::vector<double> Mass(N_Fluid_Comp,0.0);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,0,)
    {
        if(not Obstacle(i,j,k))
        for (int n = 0; n < N_Fluid_Comp; ++n)
        {
            double locMass = 0;
            for (int ii = -1; ii < 2; ii++)
            for (int jj = -1; jj < 2; jj++)
            for (int kk = -1; kk < 2; kk++)
            {
                locMass += lbPopulations(i,j,k)({n})(ii,jj,kk);
            }
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                Mass[n] += locMass;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    return Mass;
}

void FlowSolverLBM::CalcateDensityAndMomentum(void)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,lbDensity.Bcells(),)
    for (int n = 0; n < N_Fluid_Comp; ++n)
    {
        lbDensity(i,j,k)({n}) = 0.0;
        lbMomentum(i,j,k)({n}).set_to_zero();
        for (int ii = -1; ii < 2; ii++)
        for (int jj = -1; jj < 2; jj++)
        for (int kk = -1; kk < 2; kk++)
        {
            double rr = lbPopulations(i,j,k)({n})(ii,jj,kk);
            lbDensity(i,j,k)({n}) += rr;
            lbMomentum(i,j,k)({n})[0] += rr*ii;
            lbMomentum(i,j,k)({n})[1] += rr*jj;
            lbMomentum(i,j,k)({n})[2] += rr*kk;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void FlowSolverLBM::CalculateFluidVelocities(const PhaseField& Phase,
        Velocities& Vel)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,lbDensity,0,)
    if(!Obstacle(i,j,k))
    {
        const dVector3 lbFluidVelocity = (lbMomentum(i,j,k)({0}) +
                lbMomentumDelta(i,j,k)({0}) * 0.5) / lbDensity(i,j,k)({0});

        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        if(Phase.FieldsStatistics[alpha->index].State != Solid)
        {
            Vel.Phase(i,j,k)({Phase.FieldsStatistics[alpha->index].Phase}) =
                lbFluidVelocity * dxt;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Vel.CalculateAverage(Phase);
}

void FlowSolverLBM::Remesh(const int newNx, const int newNy, const int newNz,
        const BoundaryConditions& BC)
{
    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    lbDensity.Reallocate(Nx, Ny, Nz);
    Obstacle.Reallocate(Nx, Ny, Nz);
    lbMomentum.Reallocate(Nx, Ny, Nz);
    lbMomentumDelta.Reallocate(Nx, Ny, Nz);
    lbForceDensity.Reallocate(Nx, Ny, Nz);
    lbPopulations.Reallocate(Nx, Ny, Nz);
    lbPopulationsTMP.Reallocate(Nx, Ny, Nz);
}
}
