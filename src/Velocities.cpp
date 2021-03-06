#include "Velocities.h"
#include "Info.h"
#include "Settings.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "Tools/UserInterface.h"
#include "VTK.h"


namespace opensim
{
using namespace std;
/*************************************************************************/

Velocities::Velocities(const Settings& locSettings, const unsigned int boundary)
{
    this->Initialize(locSettings, boundary);
}

void Velocities::Initialize(const Settings& locSettings,
        const unsigned int boundary)
{
    thisclassname = "Velocities";

    Nx      = locSettings.Nx;
    Ny      = locSettings.Ny;
    Nz      = locSettings.Nz;
    Nphases = locSettings.Nphases;
    dx      = locSettings.dx;

    Phase.Allocate(Nx, Ny, Nz, {Nphases}, boundary);
    Average.Allocate(Nx, Ny, Nz, boundary);
    Clear();
    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void Velocities::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetXVector(Average);
    BC.SetYVector(Average);
    BC.SetZVector(Average);

    BC.SetXVector(Phase);
    BC.SetYVector(Phase);
    BC.SetZVector(Phase);
}

void Velocities::Clear()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Average,Average.Bcells(),)
    {
        Average(i,j,k).set_to_zero();
        for(int alpha = 0; alpha < Nphases; ++alpha)
        {
            Phase(i, j, k)({alpha}).set_to_zero();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void Velocities::SetAverage(dVector3& value)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Average,Average.Bcells(),)
    {
        Average(i,j,k) = value;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Velocities::SetAllPhases(dVector3& value)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase,Phase.Bcells(),)
    for(int alpha = 0; alpha < Nphases; ++alpha)
    {
        Phase(i, j, k)({alpha}) = value;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Velocities::PrescribePhaseVelocities(PhaseField& Phi)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase,Phase.Bcells(),)
    for(int alpha = 0; alpha < Nphases; ++alpha)
    {
        Phase(i, j, k)({alpha}) = Average(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Velocities::CalculateAverage(const PhaseField& Phi)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Average,Average.Bcells(),)
    {
        Average(i,j,k).set_to_zero();

        if(Phi.Interface(i,j,k))
        {
            for(auto alpha = Phi.Fields(i,j,k).cbegin();
                     alpha != Phi.Fields(i,j,k).cend(); ++alpha)
            {
                int locPindex = Phi.FieldsStatistics[alpha->index].Phase;
                Average(i,j,k) += Phase(i, j, k)({locPindex})*alpha->value;
            }
        }
        else
        {
            int pInd = Phi.FieldsStatistics[Phi.Fields(i,j,k).front().index].Phase;
            Average(i,j,k) = Phase(i,j,k)({pInd});
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
/*
void Velocities::WriteAverageVTK(int tStep)
{
    stringstream outbufer;

    outbufer << "# vtk DataFile Version 3.0\n";
    outbufer << "Velocities" << "\n";
    outbufer << "ASCII\n";
    outbufer << "DATASET RECTILINEAR_GRID\n";
    outbufer << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    outbufer << "X_COORDINATES " << Nx << " double\n";
    for (int i = 0; i < Nx; i++) outbufer << i << " ";
    outbufer << "\n";
    outbufer << "Y_COORDINATES " << Ny << " double\n";
    for (int j = 0; j < Ny; j++) outbufer << j << " ";
    outbufer << "\n";
    outbufer << "Z_COORDINATES " << Nz << " double\n";
    for (int k = 0; k < Nz; k++) outbufer << k << " ";
    outbufer << "\n";
    outbufer << "POINT_DATA " << Nx*Ny*Nz << "\n";

    for (int dir = 0; dir < 3; ++dir)
    {
        outbufer << "SCALARS v_" << dir << " double 1\n";
        outbufer << "LOOKUP_TABLE default\n";

        for (int k = 0; k < Nz; k++)
        for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++)
        {
            outbufer << Average(i,j,k)[dir] <<"\n";
        }
    }
    outbufer << "VECTORS V double\n";

    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        outbufer << Average(i,j,k)[0] << " " << Average(i,j,k)[1] << " " << Average(i,j,k)[2] << "\n";
    }

    string FileName = UserInterface::MakeFileName(VTKDir,"Velocities_", tStep, ".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbufer.rdbuf();
    vtk_file.close();
}  //  WriteVTK
*/
void Velocities::WriteAverageVTK(int tStep) const
{
    stringstream buffer;
    std::vector<int> DataTypes {PDScalars, PDVectors};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        //VTK::WriteVector(buffer, Average, "Average");
        VTK::WriteVectorAsVector(buffer, Average, "AverageV");
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "Velocities", tStep);
}  //  WriteVTK

void Velocities::WritePhaseVTK(int tStep, int idx) const
{
    stringstream buffer;
    std::vector<int> DataTypes {PDScalars, PDVectors};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        //VTK::WriteVector(buffer, Phase, "v(" + to_string(idx) + ")", idx);
        VTK::WriteVectorAsVector(buffer, Phase, "Phase", idx);
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "PhaseVelocities" + to_string(idx), tStep);
}  //  WritePhaseVTK
/*
void Velocities::WritePhaseVTK(int tStep, int idx)
{
    stringstream outbufer;

    outbufer << "# vtk DataFile Version 3.0\n";
    outbufer << "Velocities" << "\n";
    outbufer << "ASCII\n";
    outbufer << "DATASET RECTILINEAR_GRID\n";
    outbufer << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    outbufer << "X_COORDINATES " << Nx << " double\n";
    for (int i = 0; i < Nx; i++) outbufer << i << " ";
    outbufer << "\n";
    outbufer << "Y_COORDINATES " << Ny << " double\n";
    for (int j = 0; j < Ny; j++) outbufer << j << " ";
    outbufer << "\n";
    outbufer << "Z_COORDINATES " << Nz << " double\n";
    for (int k = 0; k < Nz; k++) outbufer << k << " ";
    outbufer << "\n";
    outbufer << "POINT_DATA " << Nx*Ny*Nz << "\n";

    for (int dir = 0; dir < 3; ++dir)
    {
        outbufer << "SCALARS v" << idx << "_" << dir << " double 1\n";
        outbufer << "LOOKUP_TABLE default\n";

        for (int k = 0; k < Nz; k++)
        for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++)
        {
            outbufer << Phase(i,j,k)({idx})[dir] <<"\n";
        }
    }
    outbufer << "VECTORS V" << idx << " double\n";

    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        outbufer << Phase(i,j,k)({idx})[0] << " " << Phase(i,j,k)({idx})[1] << " " << Phase(i,j,k)({idx})[2] << "\n";
    }
    stringstream converter;
    converter << "PhaseVelocities_" << idx << "_";
    string FileName = UserInterface::MakeFileName(VTKDir, converter.str(), tStep, ".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbufer.rdbuf();
    vtk_file.close();
}  //  WritePhaseVTK
*/
void Velocities::Remesh(int newNx, int newNy, int newNz, BoundaryConditions& BC)
{
    Phase.Remesh(newNx, newNy, newNz);
    Average.Remesh(newNx, newNy, newNz);

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    SetBoundaryConditions(BC);
    Info::WriteStandard(thisclassname, "Remeshed");
}

void Velocities::MoveFrame(int dx, int dy, int dz, BoundaryConditions& BC)
{
    int xBeg = (dx >= 0) + (dx < 0)*(Nx) - 1;
    int xEnd = (dx >= 0)*(Nx) + (dx < 0) - 1;
    int xInc = 1 - 2*(dx < 0);

    int yBeg = (dy >= 0) + (dy < 0)*(Ny) - 1;
    int yEnd = (dy >= 0)*(Ny) + (dy < 0) - 1;
    int yInc = 1 - 2*(dy < 0);

    int zBeg = (dz >= 0) + (dz < 0)*(Nz) - 1;
    int zEnd = (dz >= 0)*(Nz) + (dz < 0) - 1;
    int zInc = 1 - 2*(dz < 0);

    for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
    {
        Phase(i, j, k) = Phase(i + dx, j + dy, k + dz);
        Average(i, j, k) = Average(i + dx, j + dy, k + dz);
    }

    SetBoundaryConditions(BC);

    Info::WriteStandard(thisclassname, "Frame moved");

}

void Velocities::PrintPointStatistics(int x, int y, int z)
{
    stringstream message;
    message << "Point: " << x << " " << y << " " << z << endl;
    message << "Phase      Velocities " << endl;
    for (int alpha = 0; alpha < Nphases;  ++alpha)
    {
        message << alpha << ":    (" << Phase(x,y,z)({alpha})[0] << " "
                                     << Phase(x,y,z)({alpha})[1] << " "
                                     << Phase(x,y,z)({alpha})[2] << ")" << endl;
    }
    message << "Average Velocity: (" << Average(x,y,z)[0] << " "
                                     << Average(x,y,z)[1] << " "
                                     << Average(x,y,z)[2] << ")" << endl;
    Info::WriteSimple(message.str());
}

double Velocities::getMaxVelocity()
{
    maxVelocity = 0.0;
    criticalPoint[0] = -1;
    criticalPoint[1] = -1;
    criticalPoint[2] = -1;

    for(int i = 0; i < Nx; i++)
    for(int j = 0; j < Ny; j++)
    for(int k = 0; k < Nz; k++)
    for(int dir = 0; dir < 3; dir++)
    if(fabs(Average(i,j,k)[dir]) > fabs(maxVelocity))
    {
        maxVelocity = fabs(Average(i,j,k)[dir]);
        criticalPoint[0] = i;
        criticalPoint[1] = j;
        criticalPoint[2] = k;
    }
//    cout<<criticalPoint[0]<<"|"<<criticalPoint[1]<<"|"<<criticalPoint[2]<<": "<<maxVelocity<<endl;

    return maxVelocity;
}

double Velocities::getMaxVelocityOMP()
{
    int criticalPointX = -1;
    int criticalPointY = -1;
    int criticalPointZ = -1;
    double maxVelocityTemp = 0.0;

    criticalPoint[0] = -1;
    criticalPoint[1] = -1;
    criticalPoint[2] = -1;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Average,0,shared(criticalPointX, criticalPointY, criticalPointZ, maxVelocityTemp))
    for(int dir = 0; dir < 3; dir++)
    if(fabs(Average(i,j,k)[dir]) > fabs(maxVelocityTemp))
    {
        #ifdef _OPENMP
        #pragma omp critical
        #endif
        {
            // Check if still bigger
            if(fabs(Average(i,j,k)[dir]) > fabs(maxVelocityTemp))
            {
                maxVelocityTemp = fabs(Average(i,j,k)[dir]);
                criticalPointX = i;
                criticalPointY = j;
                criticalPointZ = k;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    criticalPoint[0] = criticalPointX;
    criticalPoint[1] = criticalPointY;
    criticalPoint[2] = criticalPointZ;

    maxVelocity = maxVelocityTemp;
    return maxVelocity;
}

dVector3 Velocities::getCriticalPoint()
{
    dVector3 tmp = criticalPoint;
    return tmp;
}

void Velocities::WriteStatistics(int tStep, double dt)
{
    double sim_time = tStep*dt;

    getMaxVelocity();

    ofstream output_file;
    if (tStep == 0)
    {
        output_file.open("VelocityStatistics.txt", ios::out);
        output_file << "tStep" << "\t\t\t" << "sim_time" << "\t\t\t"
                    << "x_critical\t\t\t" << "y_critical\t\t\t" << "z_critical\t\t\t"
                    << "maximal velocity" << endl;
        output_file.close();
    }

    output_file.open("VelocityStatistics.txt", ios::app);
    output_file << tStep << "\t\t\t" << sim_time << "\t\t\t"
                << criticalPoint[0] << "\t\t\t" << criticalPoint[1]<<"\t\t\t"
                << criticalPoint[2] << "\t\t\t" << maxVelocity<<endl;
    output_file.close();
}

Velocities& Velocities::operator= (const Velocities& rhs)
{
    // protect against invalid self-assignment and copy of unitialized object
    if (this != &rhs and rhs.thisclassname == "Velocities")
    {
        thisclassname = rhs.thisclassname;
        Nx = rhs.Nx;
        Ny = rhs.Ny;
        Nz = rhs.Nz;
        dx = rhs.dx;
        Nphases = rhs.Nphases;

        if (Phase.IsNotAllocated())
        {
            Phase.Allocate(Nx, Ny, Nz, {Nphases}, rhs.Phase.Bcells());
            Average.Allocate(Nx, Ny, Nz, rhs.Average.Bcells());
        }
        else if (not Phase.IsSize(rhs.Nx, rhs.Ny, rhs.Nz))
        {
            Phase.Reallocate(Nx, Ny, Nz);
            Average.Reallocate(Nx, Ny, Nz);
        }

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase,Phase.Bcells(),)
        {
            Average(i,j,k) = rhs.Average(i,j,k);
            for(int alpha = 0; alpha < Nphases; ++alpha)
            {
                Phase(i,j,k)({alpha}) = rhs.Phase(i,j,k)({alpha});
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    return *this;
}

void Velocities::Advect(BoundaryConditions& BC, double dt, int scheme)
{
    if(AverageDot.IsNotAllocated())
    {
        AverageDot.Allocate(Nx, Ny, Nz, 1);
    }
    if(not AverageDot.IsSize(Nx, Ny, Nz))
    {
        AverageDot.Reallocate(Nx, Ny, Nz);
    }
    SetBoundaryConditions(BC);

    switch(scheme)
    {
        case Upwind:
        {
            const double dx2 = 0.5/dx;

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,AverageDot,0,)
            for (auto dir = 0; dir < 3; dir++)
            {
                AverageDot(i,j,k)[dir] = dx2 *
                     ((fabs(Average(i-1,j,k)[0]) + Average(i-1,j,k)[0])*Average(i-1,j,k)[dir] +
                      (fabs(Average(i,j-1,k)[1]) + Average(i,j-1,k)[1])*Average(i,j-1,k)[dir] +
                      (fabs(Average(i,j,k-1)[2]) + Average(i,j,k-1)[2])*Average(i,j,k-1)[dir] +
                      (fabs(Average(i+1,j,k)[0]) - Average(i+1,j,k)[0])*Average(i+1,j,k)[dir] +
                      (fabs(Average(i,j+1,k)[1]) - Average(i,j+1,k)[1])*Average(i,j+1,k)[dir] +
                      (fabs(Average(i,j,k+1)[2]) - Average(i,j,k+1)[2])*Average(i,j,k+1)[dir]) -
                      (fabs(Average(i,j,k)[0]) +
                       fabs(Average(i,j,k)[1]) +
                       fabs(Average(i,j,k)[2])) * Average(i, j, k)[dir]/dx;
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            break;
        }
        default:
        {
            Info::WriteExit("Wrong advection scheme selected",
                             thisclassname, "Velocities::Advect(BC, dt, scheme)");
            exit(13);
        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,AverageDot,0,)
    for (int dir = 0; dir < 3; ++dir)
    {
        Average(i, j, k)[dir] += AverageDot(i, j, k)[dir]*dt;
        AverageDot(i, j, k)[dir] = 0.0;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    SetBoundaryConditions(BC);
}

}// namespace opensim