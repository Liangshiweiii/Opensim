#include "Orientations.h"
#include "Info.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "Settings.h"
#include "Tools/UserInterface.h"
#include "Tools/Quaternion.h"
#include "Tools.h"
#include "Crystallography.h"
#include "Velocities.h"
#include "VTK.h"

namespace opensim
{
using namespace std;

Orientations::Orientations(Settings& locSettings)
{
    Initialize(locSettings);
}

double Orientations::getMisorientation(const dMatrix3x3 RotMatA, const dMatrix3x3 RotMatB)
{
    // Taken from Kocks, Tome and Wenk
    // 'Texture and Anisotropy', Cambride University Press 1998, p. 69, eq. 9b
    return acos(((RotMatB*RotMatA.transposed()).trace()-1.0)/2.0);
}

double Orientations::getMisorientationCubic(const dMatrix3x3 RotMatA, const dMatrix3x3 RotMatB, const Crystallography& CR)
{
    // Taken from Kocks, Tome and Wenk
    // 'Texture and Anisotropy', Cambride University Press 1998, p. 69, eq. 9b

    // Only for cubic system. Needs to be generalized!

    double misorientation = 10.0;
    for (int i = 0; i < 24; i++)
    {
        double misorientationloc = acos(((CR.SymmetriesCubic[i]*RotMatB*
                RotMatA.transposed()).trace()-1.0)/2.0);
        if (abs(misorientationloc) <= abs(misorientation))
        {
            misorientation = misorientationloc;
        }
    }
    return misorientation;
}

void Orientations::Initialize(const Settings& locSettings)
{
    thisclassname = "Orientations";
    //DefaultInputFileName = ProjectInputDir + "OrientationsInput.opi";

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    dx = locSettings.dx;

    Nphases = locSettings.Nphases;

    Quaternions.Allocate(Nx, Ny, Nz, 1);
    
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Quaternions, Quaternions.Bcells(),)
    {
        Quaternions(i,j,k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    initialized = true;
    Info::WriteStandard("Orientations", "Initialized");
}

void Orientations::SetRandomGrainOrientations(PhaseField& Phase, const int seed)
{
    // To be called before EP.SetGrainsProperties(PhaseField, Orientations)
    // Warning: Will overwrite manually changed grain orientations with random angles.

	default_random_engine generator(seed);

	uniform_real_distribution<double> distributionQ1(0, Pi);
    uniform_real_distribution<double> distributionQ2(0, Pi);
    uniform_real_distribution<double> distributionQ3(0, Pi);

    for(unsigned int alpha = 0; alpha < Phase.FieldsStatistics.size(); alpha++)
    if(Phase.FieldsStatistics[alpha].Exist)
    {
        EulerAngles tempAngle;
        tempAngle.set(distributionQ1(generator),
                      distributionQ2(generator),
                      distributionQ3(generator),
                      XYZ);
        Phase.FieldsStatistics[alpha].Orientation = tempAngle.getQuaternion();
    }
}

void Orientations::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(Quaternions);
    BC.SetY(Quaternions);
    BC.SetZ(Quaternions);
}

void Orientations::Remesh(const int newNx, const int newNy, const int newNz,    // RUN LD.RemeshRotations() before executing this function!!!
        BoundaryConditions& BC)
{
    Quaternions.Remesh(newNx, newNy, newNz);
    
    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

	SetBoundaryConditions(BC);

	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Quaternions, Quaternions.Bcells(),)
	{
		Quaternions(i, j, k).setRotationMatrix();
	}
	OMP_PARALLEL_STORAGE_LOOP_END

    Info::WriteStandard(thisclassname, "Remeshed");
}

Orientations& Orientations::operator= (const Orientations& rhs)
{
    // protect against invalid self-assignment and copy of unitialized object
    if (this != &rhs and rhs.thisclassname == "Orientations")
    {
        thisclassname = rhs.thisclassname;
        //DefaultInputFileName = rhs.DefaultInputFileName;

        Nx = rhs.Nx;
        Ny = rhs.Ny;
        Nz = rhs.Nz;

        dx = rhs.dx;

        Nphases = rhs.Nphases;

        Quaternions = rhs.Quaternions;
    }
    return *this;
}

void Orientations::Write(const int tStep) const
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"Rotations_", tStep, ".dat");

    fstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be created", thisclassname, "WriteRotations");
        exit(1);
    };

    for(long int i = 0; i < Quaternions.sizeX(); ++i)
    for(long int j = 0; j < Quaternions.sizeY(); ++j)
    for(long int k = 0; k < Quaternions.sizeZ(); ++k)
    {
        for(int n = 0; n < 4; n++)
        {
            double tmp = Quaternions(i,j,k)[n];
            out.write(reinterpret_cast<char*>(&tmp), sizeof(double));
        }
    }
}

void Orientations::Read(const int tStep)
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"Rotations_",
                                                  tStep, ".dat");

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be created",
                        thisclassname, "ReadRotations");
        exit(1);
    };

    for(long int i = 0; i < Quaternions.sizeX(); ++i)
    for(long int j = 0; j < Quaternions.sizeY(); ++j)
    for(long int k = 0; k < Quaternions.sizeZ(); ++k)
    {
        double tmp[4];
            
        for(int n = 0; n < 4; n++)
        {
            inp.read(reinterpret_cast<char*>(&tmp[n]), sizeof(double));
        }        
        Quaternions(i,j,k).set(tmp[0],tmp[1],tmp[2],tmp[3]);            
    }
}

void Orientations::WriteVTK(const int tStep)
{
    stringstream buffer;
    std::vector<int> DataTypes {PDScalars};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);

    VTK::WriteQuaternion(buffer, Quaternions, "RQ");
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "Quaternions", tStep);
}

void Orientations::WriteTotalVTK(PhaseField& Phase, const int tStep) const
{
    Storage3D<dMatrix3x3, 0>    TotalRotations;
    TotalRotations.Allocate(Nx, Ny, Nz, 0);
    
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,TotalRotations,0,)
    {
        TotalRotations(i,j,k) = getTotalRotation(Phase, i, j, k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    stringstream buffer;
    std::vector<int> DataTypes {PDScalars};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        VTK::WriteMatrix(buffer, TotalRotations, "R");
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "TotalRotations", tStep);
}

void Orientations::WriteMisorientationsVTK(const int tStep, const std::string measure) const
{
    Storage3D<double, 0> tempMis;
    if(tempMis.IsNotAllocated())
    {
        tempMis.Allocate(Nx, Ny, Nz, 0);
    }
    dMatrix3x3 Unity; Unity.set_to_unity();

    double mult = 1.0;
    if (!measure.compare("deg")) mult*=180.0/Pi;

    // Calculate misorientations
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,tempMis,0,)
    {
        tempMis(i,j,k) = getMisorientation(Unity, Quaternions(i,j,k).RotationMatrix)*mult;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    
    stringstream buffer;
    std::vector<int> DataTypes {PDScalars};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        VTK::WriteScalar(buffer, tempMis, "omega");
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "Misorientations", tStep);
}

void Orientations::WriteGrainEBSDDataQuaternions(const PhaseField& Phase,
                                                        const int tStep)
{
    stringstream outbuffer;
    outbuffer << "Index\t"
              << "Phase\t"
              << "X\t"
              << "Y\t"
              << "Z\t"
              << "Quat real\t"
              << "Quat i\t"
              << "Quat j\t"
              << "Quat k\t" << endl;

    int index = 1;
    STORAGE_LOOP_BEGIN(i,j,k,Quaternions,0)
    {
        EulerAngles tempAng;
        dMatrix3x3  tempMatrix;
        Quaternion  tempQuat;

        outbuffer << index << "\t";
        if (Phase.Interface(i,j,k))
        {
            // linear interpolation in interface via quaternions
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                tempQuat += (Quaternions(i,j,k) + Phase.FieldsStatistics[alpha->index].Orientation)*alpha->value;
            }
            tempQuat.normalize();
            // Set phase index to 0 in interface
            outbuffer << 0 << "\t";
        }
        else
        {
            tempQuat = Quaternions(i,j,k) + Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Orientation;
            outbuffer << Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase + 1 << "\t";
        }

        outbuffer << i*Phase.dx << "\t"
                  << j*Phase.dx << "\t"
                  << k*Phase.dx << "\t"
                  << tempQuat[0] << "\t"
                  << tempQuat[1] << "\t"
                  << tempQuat[2] << "\t"
                  << tempQuat[3] << endl;
        index++;
    }
    STORAGE_LOOP_END
    
    string FileName = UserInterface::MakeFileName(RawDataDir, "ebsd_", tStep, ".dat");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbuffer.rdbuf();
    vtk_file.close();

    Info::WriteStandard("EBSD file", FileName);
}

void Orientations::WriteRotated100Vector(const int tStep)
{
    stringstream outbufer;
    dVector3 X;
    X[0] = 1.0;
    X[1] = 0.0;
    X[2] = 0.0;

    outbufer << "# vtk DataFile Version 3.0\n";
    outbufer << "Rotated100Vector" << "\n";
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
        outbufer << "SCALARS X_" << dir << " double 1\n";
        outbufer << "LOOKUP_TABLE default\n";

        for (int k = 0; k < Nz; k++)
        for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++)
        {
            outbufer << (Quaternions(i,j,k).RotationMatrix*X)[dir] <<"\n";
        }
    }
    outbufer << "VECTORS X double\n";

    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        outbufer << (Quaternions(i,j,k).RotationMatrix*X)[0] << " " 
                 << (Quaternions(i,j,k).RotationMatrix*X)[1] << " " 
                 << (Quaternions(i,j,k).RotationMatrix*X)[2] << "\n";
    }

    string FileName = UserInterface::MakeFileName(VTKDir,"Rotated100Vector_", tStep, ".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbufer.rdbuf();
    vtk_file.close();
}

void Orientations::Advect(Velocities& Vel, BoundaryConditions& BC, 
                                                 double dt, int scheme)
{
    if(QuaternionsDot.IsNotAllocated())
    {
        QuaternionsDot.Allocate(Nx, Ny, Nz, 1);
    }
    if(not QuaternionsDot.IsSize(Nx, Ny, Nz))
    {
        QuaternionsDot.Reallocate(Nx, Ny, Nz);
    }
    
    SetBoundaryConditions(BC);
    switch(scheme)
    {
        case Upwind:
        {                      
            const double dx  = 1.0/Vel.dx;
            const double dx2 = 0.5/Vel.dx;

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,QuaternionsDot,0,)
            {
                Quaternion val = Quaternions(i,j,k);

                Quaternion valXm = Quaternions(i-1,j,k);
                Quaternion valXp = Quaternions(i+1,j,k);
                Quaternion valYm = Quaternions(i,j-1,k);
                Quaternion valYp = Quaternions(i,j+1,k);
                Quaternion valZm = Quaternions(i,j,k-1);
                Quaternion valZp = Quaternions(i,j,k+1);

                double ux = Vel.Average(i,j,k)[0];
                double uy = Vel.Average(i,j,k)[1];
                double uz = Vel.Average(i,j,k)[2];

                double uxm = Vel.Average(i-1,j,k)[0];
                double uxp = Vel.Average(i+1,j,k)[0];

                double uym = Vel.Average(i,j-1,k)[1];
                double uyp = Vel.Average(i,j+1,k)[1];

                double uzm = Vel.Average(i,j,k-1)[2];
                double uzp = Vel.Average(i,j,k+1)[2];

                // 1st order UpWind
                QuaternionsDot(i,j,k) =
                        (valXm*(fabs(uxm) + uxm) +
                         valYm*(fabs(uym) + uym) +
                         valZm*(fabs(uzm) + uzm) +
                         valXp*(fabs(uxp) - uxp) +
                         valYp*(fabs(uyp) - uyp) +
                         valZp*(fabs(uzp) - uzp)
                        )*dx2 -
                        val*(fabs(ux) + fabs(uy) + fabs(uz))/dx;
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            break;
        }
        default:
        {
            Info::WriteExit("Wrong advection scheme selected",
                             thisclassname, "Rotations::Advect(Vel, dt, scheme)");
            exit(13);
        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,QuaternionsDot,0,)
    {
        Quaternions(i,j,k) += QuaternionsDot(i,j,k)*dt;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    
    SetBoundaryConditions(BC);
}

void Orientations::PrintPointStatistics(int x, int y, int z)
{
    cout << "Point:      (" << x << ", " << y << ", " << z << ")" << endl;
    cout << "Quaternions:\n" << Quaternions(x,y,z).print() << endl;

    cout << endl;
}

}// namespace opensim