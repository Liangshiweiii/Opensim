#define NoVoro

#include "Initializations.h"
#include "Settings.h"
#include "BoundaryConditions.h"
#include "PhaseField.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "GrainInfo.h"
#include "Info.h"
#include "Tools.h"
#include "Orientations.h"
#include "Tools/Includes.h"
#include "InterfaceEnergy.h"
#include "InterfaceMobility.h"
#include "DoubleObstacle.h"
#include "UserDrivingForce.h"
#include "DrivingForce.h"

#include <algorithm>    // std::find
#include <vector>       // std::vector

#ifndef NoVoro
#include "../tools/external/voro++/src/voro++.hh"
#endif

namespace opensim{
    double rnd(){
        return double(rand())/RAND_MAX;
    }
    using namespace std;
    #ifndef NoVoro
    using namespace voro;
    #endif

    void Initializations::CreateMicrostructure(std::string InputFileName, PhaseField& Phase, BoundaryConditions& BC, Settings& locSettings){
        string fileDir = InputFileName;
	Info::WriteLineInsert("Initialization");
	Info::WriteStandard("Source", fileDir);

	fstream inp(fileDir.c_str(), ios::in);

	if (!inp)
	{
	    Info::WriteExit("File " + fileDir + " could not be opened. Initializations::CreateMicrostructure()");
		exit(1);
	};

	int moduleLocation = UserInterface::FindModuleLocation(inp, "Initialization");
	moduleLocation = UserInterface::FindParameter(inp, moduleLocation, "InitType");
	Info::WriteLine();
	string InitType = UserInterface::ReadParameterS(inp, moduleLocation, "InitType", false);
	bool found = true;
	while(found)
	{
		/*bool*/ found = false;
		if (InitType == "Single")
		{
			int PhaseIndex = UserInterface::ReadParameterI(inp, moduleLocation, "PhaseIndex");

			Single(Phase, PhaseIndex, BC, locSettings);
			found = true;
		}

		if (InitType == "Sphere")
		{
			int PhaseIndex = UserInterface::ReadParameterI(inp, moduleLocation, "PhaseIndex");
			double Radius = UserInterface::ReadParameterD(inp, moduleLocation, "Radius");
			double x0 = UserInterface::ReadParameterD(inp, moduleLocation, "x0");
			double y0 = UserInterface::ReadParameterD(inp, moduleLocation, "y0");
			double z0 = UserInterface::ReadParameterD(inp, moduleLocation, "z0");

			Sphere(Phase, PhaseIndex, Radius, x0, y0, z0, BC, locSettings);
			found = true;
		}

		if (InitType == "Fractional")
		{
			int MajorityIndex = UserInterface::ReadParameterI(inp, moduleLocation, "MajorityIndex");
			int MinorityIndex = UserInterface::ReadParameterI(inp, moduleLocation, "MinorityIndex");
			int LayerThickness = UserInterface::ReadParameterI(inp, moduleLocation, "LayerThickness");

			Fractional(Phase, MajorityIndex, MinorityIndex, LayerThickness, BC, locSettings);
			found = true;
		}

		if (InitType == "RandomNuclei")
		{
			int phaseIndex = UserInterface::ReadParameterI(inp, moduleLocation, "PhaseIndex");
			int NParticles = UserInterface::ReadParameterI(inp, moduleLocation, "NParticles");
			int Seed = UserInterface::ReadParameterI(inp, moduleLocation, "Seed", false, -1);
			RandomNuclei(Phase, locSettings, phaseIndex, NParticles, Seed);
			found = true;
		}

		if (InitType == "QuasiRandomNuclei")
		{
			int phaseIndex = UserInterface::ReadParameterI(inp, moduleLocation, "PhaseIndex");
			int Dist = UserInterface::ReadParameterI(inp, moduleLocation, "Distance");
			int Seed = UserInterface::ReadParameterI(inp, moduleLocation, "Seed", false, -1);
			int Offset = UserInterface::ReadParameterI(inp, moduleLocation, "Offset", false, 2);
			QuasiRandomNuclei(Phase, locSettings, phaseIndex, Dist, Seed, Offset);
			found = true;
		}

		if (InitType == "Ellipsoid")
		{
			int PhaseIndex = UserInterface::ReadParameterI(inp, moduleLocation, "PhaseIndex");
			double RadiusX = UserInterface::ReadParameterD(inp, moduleLocation, "RadiusX");
			double RadiusY = UserInterface::ReadParameterD(inp, moduleLocation, "RadiusY");
			double RadiusZ = UserInterface::ReadParameterD(inp, moduleLocation, "RadiusZ");
			double x0 = UserInterface::ReadParameterD(inp, moduleLocation, "x0");
			double y0 = UserInterface::ReadParameterD(inp, moduleLocation, "y0");
			double z0 = UserInterface::ReadParameterD(inp, moduleLocation, "z0");

			Ellipsoid(Phase, PhaseIndex, RadiusX, RadiusY, RadiusZ, x0, y0, z0, BC, locSettings);
			found = true;
		}

		if (InitType == "Fractional")
		{
			int MajorityIndex = UserInterface::ReadParameterI(inp, moduleLocation, "MajorityIndex");
			int MinorityIndex = UserInterface::ReadParameterI(inp, moduleLocation, "MinorityIndex");
			double MinorityThickness = UserInterface::ReadParameterD(inp, moduleLocation, "MinorityThickness");

			Fractional(Phase, MajorityIndex, MinorityIndex, MinorityThickness, BC, locSettings);
		}

		if (InitType == "ThreeFractionals")
		{
			int Phase1Index = UserInterface::ReadParameterI(inp, moduleLocation, "Phase1Index");
			double Phase1Thickness = UserInterface::ReadParameterD(inp, moduleLocation, "Phase1Thickness");
			int Phase2Index = UserInterface::ReadParameterI(inp, moduleLocation, "Phase2Index");
			double Phase2Thickness = UserInterface::ReadParameterD(inp, moduleLocation, "Phase2Thickness");
			int Phase3Index = UserInterface::ReadParameterI(inp, moduleLocation, "Phase3Index");

			ThreeFractionals(Phase, Phase1Index, Phase1Thickness, Phase2Index, Phase2Thickness, Phase3Index, BC, locSettings);
		}

		moduleLocation = UserInterface::FindParameter(inp, moduleLocation, "InitType");
		if (moduleLocation == -1) break;

		Info::WriteLine();
		InitType = UserInterface::ReadParameterS(inp, moduleLocation, "InitType", false);
	    }
    }

    vector<vector<int> > Initializations::QuasiRandomNuclei(PhaseField& Phase, Settings& locSettings, int phaseIndex, int dist, int seed, int offset){
        vector<vector<int> > result;

	    if (seed == -1){
		srand(time(NULL));
	    }
	    else
	    {
		    srand(seed);
	    }
	
        int distx = (dist < locSettings.Nx) ? dist : 0;
        int disty = (dist < locSettings.Ny) ? dist : 0;
        int distz = (dist < locSettings.Nz) ? dist : 0;
        double threshold = 0.05;

        for (int i = distx; i < locSettings.Nx; i += 2 * dist)
        for (int j = disty; j < locSettings.Ny; j += 2 * dist)
        for (int k = distz; k < locSettings.Nz; k += 2 * dist)
        {
            double chance = double(rand()) / double(RAND_MAX);
            if (chance > threshold)
            {
                int di = (distx == 0) ? 0 : (rand() % (2 * offset) - offset);
                int dj = (disty == 0) ? 0 : (rand() % (2 * offset) - offset);
                int dk = (distz == 0) ? 0 : (rand() % (2 * offset) - offset);
                if (i + di >= 0 and i + di < locSettings.Nx and
                    j + dj >= 0 and j + dj < locSettings.Ny and
                    k + dk >= 0 and k + dk < locSettings.Nz)
                {
                    Phase.PlantGrainNucleus(phaseIndex, i + di, j + dj, k + dk);

                    vector<int> temp(3, 0);
                    temp[0] = i + di;
                    temp[1] = j + dj;
                    temp[2] = k + dk;

                    result.push_back(temp);
                }
            }
        }

        return result;
        }

        void Initializations::RandomNuclei(PhaseField& Phase, Settings& locSettings, int phaseIndex, int Nparticles, int seed)
        {
        if (seed == -1)
        {
            srand(time(NULL));
        }
        else
        {
            srand(seed);
        }
        for (int n = 0; n < Nparticles; n++)
        {
            bool planted = false;
            int iterations = 0;
            while(!planted and iterations < 1000)
            {
                iterations++;
                int di = rand() % locSettings.Nx;
                int dj = rand() % locSettings.Ny;
                int dk = rand() % locSettings.Nz;
                
                bool freecell = true;

                for (int ii = - locSettings.iWidth; ii <= locSettings.iWidth; ii++)
                for (int jj = - locSettings.iWidth; jj <= locSettings.iWidth; jj++)
                for (int kk = - locSettings.iWidth; kk <= locSettings.iWidth; kk++)
                if ((di+ii > -Phase.Fields.Bcells() and di+ii < locSettings.Nx + Phase.Fields.Bcells() and
                    dj+jj > -Phase.Fields.Bcells() and dj+jj < locSettings.Ny + Phase.Fields.Bcells() and
                    dk+kk > -Phase.Fields.Bcells() and dk+kk < locSettings.Nz + Phase.Fields.Bcells()) and		
                    Phase.Fields(di+ii,dj+jj,dk+kk).flag)
                {
                    freecell = false;
                    break;
                }
                if(freecell)
                {
                    int locIndex = Phase.PlantGrainNucleus(phaseIndex, di, dj, dk);
                    int Q1 = rand() % 180;
                    int Q2 = rand() % 180;
                    int Q3 = rand() % 180;
                    EulerAngles locAngles({Q1*Pi/180, Q2*Pi/180, Q3*Pi/180}, XYZ);
                    Phase.FieldsStatistics[locIndex].Orientation = locAngles.getQuaternion();
                    planted = true;
                }
            }
        }
    }
    /**
     * 在平面的反面初始化一个新的相Initializes a new phase on the negative side of a plane
     * @param Phase       Phase field相场
     * @param Point       Point on plane平面上的点
     * @param Orientation Orientation of the plane平面方向
     * @param PhaseIndex  Phase index序参数
     * @param BC          Boundary conditions边界条件
     * @param locSettings Project settings
     */
    int Initializations::SectionalPlane(PhaseField& Phase, const int PhaseIndex,
        const dVector3 Point, const dVector3 Orientation,
        const BoundaryConditions& BC, const Settings& locSettings, const bool NewGrain, const int FieldIndex){

        const dVector3 n        = Orientation.normalized();
        const double   iWidth   = locSettings.iWidth;

        int locIndex;
        if (NewGrain)
        {
            locIndex = Phase.AddGrainInfo(PhaseIndex);
        }
        else
        {
            locIndex = FieldIndex;
        }

        STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
        {
            const dVector3 x = {double(i),double(j),double(k)}; // 空间位置position in space
            const double distance = (x-Point)*n;   // 离平面的距离distance from plane

            // 如果点在反面if the point is on the negative side
            if (distance <= 0.5*iWidth)
            {
                // 如果点不在界面上if the point is not in the interface
                if (distance < - 0.5*iWidth)
                {
                    Phase.Fields(i, j, k).clear();
                    Phase.Fields(i, j, k).set(locIndex, 1.0);
                    Phase.Fields(i,j,k).flag = 0;
                }
                else  // 如果点在界面上if the point is in the interface
                {
                    const double IntProf = 0.5 - 0.5*sin(Pi*(distance)/iWidth);
                    for(auto alpha = Phase.Fields(i,j,k).begin();
                            alpha < Phase.Fields(i,j,k).end(); alpha++)
                    {
                        alpha->value *= 1.0 - IntProf;
                    }
                    Phase.Fields(i, j, k).set(locIndex, IntProf);
                    Phase.Fields(i,j,k).flag = 2;
                }
            }
        }
        STORAGE_LOOP_END

        Phase.Finalize(BC);
        return locIndex;
    }

    int Initializations::Sphere(PhaseField& Phase, int PhaseIndex,
                             double Radius, double x0, double y0, double z0,
                             BoundaryConditions& BC, Settings& locSettings)
{
    double iWidth = locSettings.iWidth;
    int locIndex = Phase.AddGrainInfo(PhaseIndex);

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        double rad = sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0)+(k-z0)*(k-z0));
        if (rad < Radius - iWidth*0.5)
        {
            Phase.Fields(i, j, k).clear();
            Phase.Fields(i, j, k).set(locIndex, 1.0);
            Phase.Fields(i,j,k).flag = 0;
        }
        else if (rad < Radius + iWidth*0.5)
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(rad - Radius)/iWidth);
            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha < Phase.Fields(i,j,k).end(); alpha++)
            {
                alpha->value *= 1.0 - IntProf;
            }
            Phase.Fields(i, j, k).set(locIndex, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.Finalize(BC);
    return locIndex;
}
    void Initializations::Young4Periodic(PhaseField& Phase, int PhaseIndex,
        BoundaryConditions& BC, Settings& locSettings)
{
    int Nx = locSettings.Nx;
    int Ny = locSettings.Ny;
    int Nz = locSettings.Nz;

    int locIndex0 = Phase.AddGrainInfo(PhaseIndex);
    int locIndex1 = Phase.AddGrainInfo(PhaseIndex+1);
    int locIndex2 = Phase.AddGrainInfo(PhaseIndex+2);
    int locIndex3 = Phase.AddGrainInfo(PhaseIndex+3);
    int locIndex4 = Phase.AddGrainInfo(PhaseIndex+4);
    int locIndex5 = Phase.AddGrainInfo(PhaseIndex+5);
    int locIndex6 = Phase.AddGrainInfo(PhaseIndex+6);
    int locIndex7 = Phase.AddGrainInfo(PhaseIndex+7);

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i,j,k).flag = 2;

        if (j <= Ny/2)
        {
            if      (i <= Nx/2 and k >  Nz/4 and k <= Nz*3./4)
                Phase.Fields(i, j, k).set(locIndex1,1);
            else if (i >  Nx/2 and k >  Nz/4 and k <= Nz*3./4)
                Phase.Fields(i, j, k).set(locIndex2,1);
            else if (i >  Nx/4 and i <= Nx*3./4 and k <= Nz/4)
                Phase.Fields(i, j, k).set(locIndex3,1);
            else if (i >  Nx/4 and i <= Nx*3./4 and k >  Nz*3./4)
                Phase.Fields(i, j, k).set(locIndex3,1);
            else Phase.Fields(i,j,k).set(locIndex0,1);
        }
        else
        {
            if      (i <= Nx/8 and k >  Nz/8 and k <= Nz*5./8)
                Phase.Fields(i, j, k).set(locIndex5,1);
            else if (i >  Nx/8 and i <= Nx*5./8 and k >  Nz/8 and k <= Nz*5./8)
                Phase.Fields(i, j, k).set(locIndex4,1);
            else if (i >  Nx*5./8 and k >  Nz/8 and k <= Nz*5./8)
                Phase.Fields(i, j, k).set(locIndex5,1);
            else if (i >  Nx*3./8 and i <= Nx*7./8 and k <= Nz/8)
                Phase.Fields(i, j, k).set(locIndex6,1);
            else if (i >  Nx*3./8 and i <= Nx*7./8 and k >  Nz*5./8)
                Phase.Fields(i, j, k).set(locIndex6,1);
            else Phase.Fields(i,j,k).set(locIndex7,1);
        }
    }
    STORAGE_LOOP_END

    Phase.SetBoundaryConditions(BC);
    Phase.CalculateLaplacians();
}

void Initializations::Young4(PhaseField& Phase, int PhaseIndex,
        BoundaryConditions& BC, Settings& locSettings)
{
    int Nx = locSettings.Nx;
    int Ny = locSettings.Ny;
    int Nz = locSettings.Nz;

    int locIndex0 = Phase.AddGrainInfo(PhaseIndex);
    int locIndex1 = Phase.AddGrainInfo(PhaseIndex);
    int locIndex2 = Phase.AddGrainInfo(PhaseIndex);
    int locIndex3 = Phase.AddGrainInfo(PhaseIndex);

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i,j,k).flag = 2;

        if(i <  (Nx/4))                           Phase.Fields(i, j, k).set(locIndex0,1);
        if(i >= (Nx/4) && j <  Ny/3)              Phase.Fields(i, j, k).set(locIndex1,1);
        if(i >= (Nx/4) && j >= Ny/3 && k <  Nz/2) Phase.Fields(i, j, k).set(locIndex2,1);
        if(i >= (Nx/4) && j >= Ny/3 && k >= Nz/2) Phase.Fields(i, j, k).set(locIndex3,1);
    }
    STORAGE_LOOP_END

    Phase.SetBoundaryConditions(BC);
    Phase.CalculateLaplacians();
}

void Initializations::Young3(PhaseField& Phase, int PhaseIndex, BoundaryConditions& BC, Settings& locSettings)
{
    int Nx = locSettings.Nx;
    int Nz = locSettings.Nz;

    int locIndex0 = Phase.AddGrainInfo(0);
    int locIndex1 = Phase.AddGrainInfo(1);
    int locIndex2 = Phase.AddGrainInfo(2);
    int locIndex3 = Phase.AddGrainInfo(3);

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i,j,k).clear();
        Phase.Fields(i,j,k).flag = 2;
        if(i < (Nx/2)) Phase.Fields(i, j, k).set(locIndex0,1);
        else           Phase.Fields(i, j, k).set(locIndex1,1);
    }
    STORAGE_LOOP_END

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    if(k < (Nz/4) || k >= (Nz*3/4))
    {
        Phase.Fields(i,j,k).clear();

        if(i <  (Nx/4) || i >= (Nx*3/4)) Phase.Fields(i, j, k).set(locIndex2,1);
        else Phase.Fields(i, j, k).set(locIndex3,1);

        Phase.Fields(i,j,k).flag = 2;
    }
    STORAGE_LOOP_END

    Phase.SetBoundaryConditions(BC);
    Phase.CalculateLaplacians();
}

int Initializations::Single(PhaseField& Phase, int PhaseIndex, BoundaryConditions& BC, Settings& locSettings)
{
    int locIndex = Phase.AddGrainInfo(PhaseIndex);
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();
        Phase.Fields(i, j, k).set(locIndex, 1.0);
    }
    STORAGE_LOOP_END

    Phase.Finalize(BC);
    return locIndex;
}

int Initializations::Zlayer(PhaseField& Phase, int PhaseIndex, int Position,
                            int LayerThickness, BoundaryConditions& BC, Settings& locSettings)
{
    double iWidth = locSettings.iWidth;

    int index = Phase.AddGrainInfo(PhaseIndex);

    int Zmin = Position - LayerThickness/2;
    int Zmax = Position + LayerThickness/2;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        if (k > Zmin + iWidth/2 and k < Zmax - iWidth/2)
        {
            Phase.Fields(i, j, k).clear();
            Phase.Fields(i, j, k).set(index, 1);
        }
        if (k >= Zmin - iWidth/2 and k <= Zmin + iWidth/2)
        {
            double IntProfile = 0.5 - 0.5*sin(Pi*(k - Zmin)/iWidth);
            for(auto beta = Phase.Fields(i,j,k).begin();
                     beta != Phase.Fields(i,j,k).end(); beta++)
            {
                beta->value *= (IntProfile);
            }
            Phase.Fields(i, j, k).set(index, 1.0 - IntProfile);
            //Phase.Fields(i,j,k).flag = 2;
        }
        if (k >= Zmax - iWidth/2 and k <= Zmax + iWidth/2)
        {
            double IntProfile = 0.5 - 0.5*sin(Pi*(k - Zmax)/iWidth);
            for(auto beta = Phase.Fields(i,j,k).begin();
                     beta != Phase.Fields(i,j,k).end(); beta++)
            {
                beta->value *= 1.0 - IntProfile;
            }
            Phase.Fields(i, j, k).set(index, IntProfile);
            //Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.Finalize(BC);
    return index;
}
//Fractions---两组分体系的初始化。
vector<int> Initializations::Fractional(PhaseField& Phase,
                                 int MajorityPhaseIndex, int MinorityPhaseIndex,
                                 double MinorityPhaseLayerThickness,
                                 BoundaryConditions& BC, Settings& locSettings)
{

    if( BC.BC0Z != NoFlux || BC.BCNZ != NoFlux )                  //判断是否有一个边界不为Noflux类型
    {
        Info::WriteWarning(
                "Consider NoFlux Boundary condition at least along \n"
                "Z-direction for Fractional initialization", "Initializations", "Fractional");
    }

    double iWidth = locSettings.iWidth;              //界面宽度初始化

    int index1 = Phase.AddGrainInfo(MajorityPhaseIndex);             //相的存储的元素数目
    int index2 = Phase.AddGrainInfo(MinorityPhaseIndex);             //相的存储的元素数目

    int offset = MinorityPhaseLayerThickness;                        //劣势相的相层厚度

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();

        if (k > offset + iWidth/2)
        {
            Phase.Fields(i, j, k).set(index1, 1);
        }
        else if (k < offset - iWidth/2)
        {
            Phase.Fields(i, j, k).set(index2, 1);
        }
        else
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(k - offset)/iWidth);
            Phase.Fields(i, j, k).set(index1, 1.0 - IntProf);
            Phase.Fields(i, j, k).set(index2, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }                   //使用宏来初始化结构体
    STORAGE_LOOP_END

    Phase.Finalize(BC);
    return vector<int>{index1,index2};
}

vector<int> Initializations::ThreeFractionals(PhaseField& Phase,
                                 int MajorityPhaseIndex,
                                 double MajorityPhaseLayerThickness,
                                 int MinorityPhaseIndex1,
                                 double MinorityPhaseLayerThickness1,
                                 int MinorityPhaseIndex2,
                                 BoundaryConditions& BC, Settings& locSettings)
{

    if( BC.BC0Z != NoFlux || BC.BCNZ != NoFlux )
    {
        Info::WriteWarning(
                "Consider NoFlux Boundary condition atleast along \n"
                "Z-direction for Fractional initialization", "Initializations", "Fractional");
    }

    double iWidth = locSettings.iWidth;

    int index1 = Phase.AddGrainInfo(MajorityPhaseIndex);
    int index2 = Phase.AddGrainInfo(MinorityPhaseIndex1);
    int index3 = Phase.AddGrainInfo(MinorityPhaseIndex2);

    int offset1 = MajorityPhaseLayerThickness;
    int offset2 = MinorityPhaseLayerThickness1;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();
        Phase.Fields(i,j,k).flag = 0;

        if (k < offset1 - iWidth/2)
        {
            Phase.Fields(i, j, k).set(index1, 1);
        }
        else if ((k > offset1 + iWidth/2) and (k < offset1+offset2 - iWidth/2))
        {
            Phase.Fields(i, j, k).set(index2, 1);
        }
        else if (k > offset1+offset2 + iWidth/2)
        {
            Phase.Fields(i, j, k).set(index3, 1);
        }
        else if ((k <= offset1 + iWidth/2) and (k >= offset1 - iWidth/2))
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(k - offset1)/iWidth);
            Phase.Fields(i, j, k).set(index2, 1.0 - IntProf);
            Phase.Fields(i, j, k).set(index1, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
        else if ((k <= offset1+offset2 + iWidth/2) and
                 (k >= offset1+offset2 - iWidth/2))
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(k - (offset1+offset2))/iWidth);
            Phase.Fields(i, j, k).set(index3, 1.0 - IntProf);
            Phase.Fields(i, j, k).set(index2, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.Finalize(BC);
    return vector<int>{index1,index2,index3};
}

vector<int> Initializations::FractionalSharpInterface(PhaseField& Phase,
                                 int MajorityPhaseIndex, int MinorityPhaseIndex,
                                 double MinorityPhaseLayerThickness,
                                 BoundaryConditions& BC, Settings& locSettings)
{
    int index1 = Phase.AddGrainInfo(MajorityPhaseIndex);
    int index2 = Phase.AddGrainInfo(MinorityPhaseIndex);

    int offset = MinorityPhaseLayerThickness;
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();

        if (i <= offset)
        {
            Phase.Fields(i, j, k).set(index2, 1);
        }
        else
        {
            Phase.Fields(i, j, k).set(index1, 1);
        }
        Phase.Fields(i,j,k).flag = 2;
    }
    STORAGE_LOOP_END

    Phase.SetBoundaryConditions(BC);
    Phase.CalculateLaplacians();
    return vector<int>{index1,index2};
}

vector<int> Initializations::TwoWalls(PhaseField& Phase, int ChannelPhaseIndex, int WallsPhaseIndex,
                                 double WallsThickness, BoundaryConditions& BC, Settings& locSettings)
{
    int Nz = locSettings.Nz;
    double iWidth = locSettings.iWidth;

    int index1 = Phase.AddGrainInfo(ChannelPhaseIndex);
    int index2 = Phase.AddGrainInfo(WallsPhaseIndex);

    int offset = WallsThickness;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();
        Phase.Fields(i,j,k).flag = 0;

        if (k > offset + iWidth/2 && k < Nz - offset - 1 - iWidth/2)
        {
            Phase.Fields(i, j, k).set(index1, 1);
        }
        else if (k < offset - iWidth/2 || k > Nz - offset - 1 + iWidth/2)
        {
            Phase.Fields(i, j, k).set(index2, 1);
        }
        else if (k >= offset - iWidth/2 && k <= offset + iWidth/2)
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(k - offset)/iWidth);
            Phase.Fields(i, j, k).set(index1, 1.0 - IntProf);
            Phase.Fields(i, j, k).set(index2, IntProf);

            Phase.Fields(i,j,k).flag = 2;
        }
        else if (k >= Nz - offset - 1 - iWidth/2 && k <= Nz - offset - 1 + iWidth/2)
        {
            double IntProf = 0.5 + 0.5*sin(Pi*(Nz - k - offset - 1)/iWidth);
            Phase.Fields(i, j, k).set(index1, IntProf);
            Phase.Fields(i, j, k).set(index2, 1.0 - IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.Finalize(BC);

    return vector<int>{index1,index2};
}

vector<int> Initializations::TwoDifferentWalls(PhaseField& Phase, int ChannelPhaseIndex, int WallsPhaseIndex,
                                 double WallsThickness, BoundaryConditions& BC, Settings& locSettings)
{
    int Nz = locSettings.Nz;
    double iWidth = locSettings.iWidth;

    int index1 = Phase.AddGrainInfo(ChannelPhaseIndex);
    int index2 = Phase.AddGrainInfo(WallsPhaseIndex);
    int index3 = Phase.AddGrainInfo(WallsPhaseIndex);
    int offset = WallsThickness;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();

        if (k > offset + iWidth/2 && k < (Nz) - offset - iWidth/2 - 1)
        {
            Phase.Fields(i, j, k).set(index1, 1);
        }
        else if (k < offset - iWidth/2)
        {
            Phase.Fields(i, j, k).set(index2, 1);
        }
        else if (k > (Nz) - offset + iWidth/2 - 1)
        {
            Phase.Fields(i, j, k).set(index3, 1);
        }
        else if (k >= offset - iWidth/2 && k <= offset + iWidth/2)
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(k - offset)/iWidth);
            Phase.Fields(i, j, k).set(index1, 1.0 - IntProf);
            Phase.Fields(i, j, k).set(index2, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
        else if (k >= (Nz) - offset - iWidth/2 - 1 && k <= (Nz) - offset + iWidth/2 - 1)
        {
            double IntProf = 0.5 + 0.5*sin(Pi*((Nz) - k - offset - 1)/iWidth);
            Phase.Fields(i, j, k).set(index1, IntProf);
            Phase.Fields(i, j, k).set(index3, 1.0 - IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.Finalize(BC);

    return vector<int>{index1,index2,index3};
}

double PinE(const double x, const double y, const double z, const double a, const double b, const double c) /// Auxiliary function, check whether a point (x,y,z) lies within ellipsoid with axes a, b, c
{
  return (x*x/a/a + y*y/b/b + z*z/c/c - 1.0);
}

int Initializations::Ellipsoid(PhaseField& Phase,
                                int PhaseIndex,
                                double RadiusX, double RadiusY, double RadiusZ,
                                double x0, double y0, double z0,
                                BoundaryConditions& BC, Settings& locSettings)
{
    double iWidth = locSettings.iWidth;
    int index = Phase.AddGrainInfo(PhaseIndex);

    ///< Inner ellipsoid
    double Rxi = RadiusX - iWidth*0.5;
    double Ryi = RadiusY - iWidth*0.5;
    double Rzi = RadiusZ - iWidth*0.5;
    ///< Outer ellipsoid
    double Rxo = RadiusX + iWidth*0.5;
    double Ryo = RadiusY + iWidth*0.5;
    double Rzo = RadiusZ + iWidth*0.5;

    double Rxi2 = Rxi*Rxi;
    double Ryi2 = Ryi*Ryi;
    double Rzi2 = Rzi*Rzi;

    double Rxo2 = Rxo*Rxo;
    double Ryo2 = Ryo*Ryo;
    double Rzo2 = Rzo*Rzo;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        double radX = (i-x0);
        double radY = (j-y0);
        double radZ = (k-z0);

        if ((radX*radX/Rxi2 + radY*radY/Ryi2 + radZ*radZ/Rzi2) <= 1.0)  ///< Point in the pure phase
        {
            Phase.Fields(i, j, k).clear();     ///< remove all phase fields that may already be present
            Phase.Fields(i, j, k).set(index, 1.0);
        }
        else if ((radX*radX/Rxo2 + radY*radY/Ryo2 + radZ*radZ/Rzo2) <= 1.0)  ///< Point in the interface
        {
            // Find coordinates inside the inteface
            double r1 = -iWidth*0.5;
            double r2 =  iWidth*0.5;
            double rr = (r1 + r2)*0.5;
            double tolerance = 1e-6;
            while (fabs(r2 - r1) > tolerance)
            {
                if(PinE(radX, radY, radZ, RadiusX+r1, RadiusY+r1, RadiusZ+r1) > 0.0)
                   r1 += 0.25*(r2 - r1);
                else r1 -= 0.25*(r2 - r1);
                if(PinE(radX, radY, radZ, RadiusX+r2, RadiusY+r2, RadiusZ+r2) < 0.0)
                    r2 -= 0.25*(r2 - r1);
                else r2 += 0.25*(r2 - r1);

                rr = (r1 + r2)*0.5;
            }

            double IntProf = 0.5 - 0.5*sin(Pi*rr/iWidth);

            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha < Phase.Fields(i,j,k).end(); ++alpha)
            {
                alpha->value *= 1.0 - IntProf;
            }
            Phase.Fields(i, j, k).set(index, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.Finalize(BC);

    return index;
}

void Initializations::Read(PhaseField& Phase, string FileName, BoundaryConditions& BC, Settings& locSettings)
{
    Phase.Read(FileName, BC);

    Phase.FieldsStatistics[0].Exist = true;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        if (Phase.Fields(i,j,k).size() == 0)
        {
            Phase.Fields(i,j,k).set(0, 1.0);
        }
        Phase.Fields(i,j,k).flag = 2;
    }
    STORAGE_LOOP_END

    Phase.SetBoundaryConditions(BC);
    Phase.CalculateLaplacians();
}

int Initializations::CubeSimple(PhaseField& Phase, const int PhaseIndex,
        const double Lx, const double Ly, const double Lz,
        const int i0, const int j0, const int k0,
        const BoundaryConditions& BC, const Settings& locSettings)
{
    const double iWidth = locSettings.iWidth;

    const int index = Phase.AddGrainInfo(PhaseIndex);
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        // Determine the area where the phase-field is one
        const int SizeX = (Lx == locSettings.Nx) ? Lx : (Lx/2 - iWidth/2);
        const int SizeY = (Ly == locSettings.Ny) ? Ly : (Ly/2 - iWidth/2);
        const int SizeZ = (Lz == locSettings.Nz) ? Lz : (Lz/2 - iWidth/2);

        Phase.Fields(i,j,k).flag = 2;
        // Determine relative coordinates
        int ii = abs(i - i0) - SizeX;
        int jj = abs(j - j0) - SizeY;
        int kk = abs(k - k0) - SizeZ;
        // Set cubic phase field
        if (ii < 0 and jj < 0 and kk < 0)
        {
            Phase.Fields(i, j, k).clear();
            Phase.Fields(i, j, k).set(index, 1);
        }
        else if (ii <= iWidth and jj <= iWidth and kk <= iWidth)
        {
            const double xx = max(ii, max(jj, kk));
            const double IntProf = 0.5 + 0.5 * cos(Pi * xx / iWidth);
            for(auto alpha = Phase.Fields(i,j,k).begin();
                     alpha < Phase.Fields(i,j,k).end(); alpha++)
            {
                alpha->value *= 1.0 - IntProf;
            }
            Phase.Fields(i,j,k).set(index, IntProf);
        }
    }
    STORAGE_LOOP_END

    Phase.Finalize(BC);

    return index;
}

void Initializations::Cube(PhaseField& Phase, int PFIndex,
                             double Lx, double Ly, double Lz,
                             double x0, double y0, double z0,
                             BoundaryConditions& BC, Settings& locSettings)
{
    int iWidth = locSettings.iWidth;
    int halfIWidth = iWidth/2;
    int Nx = locSettings.Nx;
    int Ny = locSettings.Ny;
    int Nz = locSettings.Nz;


//    int locIndex = Phase.AddGrainInfo(PhaseIndex);
//    int grainIndex = 1;
//    int locIndex = Phase.FieldsStatistics.add_phase(grainIndex,PhaseIndex);
    int locIndex = PFIndex;
    //corners
//    SphereFixedIdx(Phase, PhaseIndex, halfIWidth,
//            x0 - Lx*0.5 /*+ halfIWidth*/, y0 - Ly*0.5 /*+ halfIWidth*/, z0 - Lz*0.5 /*+ halfIWidth*/, locSettings);
//    SphereFixedIdx(Phase, PhaseIndex, halfIWidth,
//            x0 - Lx*0.5 /*+ halfIWidth*/, y0 + Ly*0.5 /*- halfIWidth*/, z0 - Lz*0.5 /*+ halfIWidth*/, locSettings);
//    SphereFixedIdx(Phase, PhaseIndex, halfIWidth,
//            x0 + Lx*0.5 /*- halfIWidth*/, y0 + Ly*0.5 /*- halfIWidth*/, z0 - Lz*0.5 /*+ halfIWidth*/, locSettings);
//    SphereFixedIdx(Phase, PhaseIndex, halfIWidth,
//            x0 + Lx*0.5 /*- halfIWidth*/, y0 - Ly*0.5 /*+ halfIWidth*/, z0 - Lz*0.5 /*+ halfIWidth*/, locSettings);
//    SphereFixedIdx(Phase, PhaseIndex, halfIWidth,
//            x0 - Lx*0.5 /*+ halfIWidth*/, y0 - Ly*0.5 /*+ halfIWidth*/, z0 + Lz*0.5 /*- halfIWidth*/, locSettings);
//    SphereFixedIdx(Phase, PhaseIndex, halfIWidth,
//            x0 - Lx*0.5 /*+ halfIWidth*/, y0 + Ly*0.5 /*- halfIWidth*/, z0 + Lz*0.5 /*- halfIWidth*/, locSettings);
//    SphereFixedIdx(Phase, PhaseIndex, halfIWidth,
//            x0 + Lx*0.5 /*- halfIWidth*/, y0 + Ly*0.5 /*- halfIWidth*/, z0 + Lz*0.5 /*- halfIWidth*/, locSettings);
//    SphereFixedIdx(Phase, PhaseIndex, halfIWidth,
//            x0 + Lx*0.5 /*- halfIWidth*/, y0 - Ly*0.5 /*+ halfIWidth*/, z0 + Lz*0.5/* - halfIWidth*/, locSettings);

    //edges
    Cylinder(Phase, locIndex,halfIWidth,Lx,0,
             x0,
             y0+Ly*0.5-halfIWidth-1,
             z0+Lz*0.5-halfIWidth-1,
             locSettings);
    Cylinder(Phase, locIndex,halfIWidth,Lx,0,
             x0,
             y0-Ly*0.5+halfIWidth+1,
             z0+Lz*0.5-halfIWidth-1,
             locSettings);
    Cylinder(Phase, locIndex,halfIWidth,Lx,0,
             x0,
             y0+Ly*0.5-halfIWidth-1,
             z0-Lz*0.5+halfIWidth+1,
             locSettings);
    Cylinder(Phase, locIndex,halfIWidth,Lx,0,
             x0,
             y0-Ly*0.5+halfIWidth+1,
             z0-Lz*0.5+halfIWidth+1,
             locSettings);

    Cylinder(Phase, locIndex,halfIWidth,Ly,1,
             x0+Lx*0.5-halfIWidth-1,
             y0,
             z0+Lz*0.5-halfIWidth-1,
             locSettings);
    Cylinder(Phase, locIndex,halfIWidth,Ly,1,
             x0-Lx*0.5+halfIWidth+1,
             y0,
             z0+Lz*0.5-halfIWidth-1,
             locSettings);
    Cylinder(Phase, locIndex,halfIWidth,Ly,1,
             x0+Lx*0.5-halfIWidth-1,
             y0,
             z0-Lz*0.5+halfIWidth+1,
             locSettings);
    Cylinder(Phase, locIndex,halfIWidth,Ly,1,
             x0-Lx*0.5+halfIWidth+1,
             y0,
             z0-Lz*0.5+halfIWidth+1,
             locSettings);

    Cylinder(Phase, locIndex,halfIWidth,Lz,2,
             x0+Lx*0.5-halfIWidth-1,
             y0+Ly*0.5-halfIWidth-1,
             z0,
             locSettings);
    Cylinder(Phase, locIndex,halfIWidth,Lz,2,
             x0-Lx*0.5+halfIWidth+1,
             y0+Ly*0.5-halfIWidth-1,
             z0,
             locSettings);
    Cylinder(Phase, locIndex,halfIWidth,Lz,2,
             x0+Lx*0.5-halfIWidth-1,
             y0-Ly*0.5+halfIWidth+1,
             z0,
             locSettings);
    Cylinder(Phase, locIndex,halfIWidth,Lz,2,
             x0-Lx*0.5+halfIWidth+1,
             y0-Ly*0.5+halfIWidth+1,
             z0,
             locSettings);

    //inner part
    for(int i = x0 - Lx/2 + halfIWidth; i <= x0 + Lx/2 - halfIWidth; ++i)
    for(int j = y0 - Ly/2 + halfIWidth; j <= y0 + Ly/2 - halfIWidth; ++j)
    for(int k = z0 - Lz/2 + halfIWidth; k <= z0 + Lz/2 - halfIWidth; ++k)
    {
//        if(fabs(i-x0) <= Lx*0.5 - halfIWidth and
//           fabs(j-y0) <= Ly*0.5 - halfIWidth and
//           fabs(k-z0) <= Lz*0.5 - halfIWidth)
        {
            Phase.Fields(i,j,k).clear();
            Phase.Fields(i,j,k).set(locIndex,1.0);
        }
    }

    //faces
    for(int n = -halfIWidth; n <= halfIWidth; n ++)
    {
        double IntProf = 0.5 - 0.5*sin(Pi*n/(iWidth));

        for(int i = 0; i < Nx; ++i)
        for(int j = 0; j < Ny; ++j)
        {
            if(fabs(i-x0) <= Lx*0.5 - halfIWidth and fabs(j-y0) <= Ly*0.5 - halfIWidth)
            {
                int k = z0 + Lz*0.5 + n;
                if(IntProf > Phase.Fields(i, j, k)[locIndex])
                {
                    for(auto alpha = Phase.Fields(i,j,k).begin();
                             alpha < Phase.Fields(i,j,k).end(); alpha++)
                    {
                        if(alpha->index != locIndex)
                            alpha->value *= 1.0 - IntProf;
                    }
                    Phase.Fields(i, j, k).set(locIndex,IntProf);
                    Phase.Fields(i,j,k).flag = 2;
                }

                k = z0 - Lz*0.5 - n;
                if(IntProf > Phase.Fields(i, j, k)[locIndex])
                {
                    for(auto alpha = Phase.Fields(i,j,k).begin();
                             alpha < Phase.Fields(i,j,k).end(); alpha++)
                    {
                        if(alpha->index != locIndex)
                            alpha->value *= 1.0 - IntProf;
                    }
                    Phase.Fields(i, j, k).set(locIndex,IntProf);
                    Phase.Fields(i,j,k).flag = 2;
                }
            }
        }

        for(int i = 0; i < Nx; ++i)
        for(int k = 0; k < Ny; ++k)
        {
            if(fabs(i-x0) <= Lx*0.5 - halfIWidth and fabs(k-z0) <= Lz*0.5 - halfIWidth)
            {
                int j = y0 + Ly*0.5 + n;
                if(IntProf > Phase.Fields(i, j, k)[locIndex])
                {
                    for(auto alpha = Phase.Fields(i,j,k).begin();
                             alpha < Phase.Fields(i,j,k).end(); alpha++)
                    {
                        if(alpha->index != locIndex)
                            alpha->value *= 1.0 - IntProf;
                    }
                    Phase.Fields(i, j, k).set(locIndex,IntProf);
                    Phase.Fields(i,j,k).flag = 2;
                }

                j = y0 - Ly*0.5 - n;
                if(IntProf > Phase.Fields(i, j, k)[locIndex])
                {
                    for(auto alpha = Phase.Fields(i,j,k).begin();
                             alpha < Phase.Fields(i,j,k).end(); alpha++)
                    {
                        if(alpha->index != locIndex)
                            alpha->value *= 1.0 - IntProf;
                    }
                    Phase.Fields(i, j, k).set(locIndex,IntProf);
                    Phase.Fields(i,j,k).flag = 2;
                }
            }
        }

        for(int j = 0; j < Ny; ++j)
        for(int k = 0; k < Nz; ++k)
        {
            if(fabs(j-y0) <= Ly*0.5 - halfIWidth and fabs(k-z0) <= Lz*0.5 - halfIWidth)
            {
                int i = x0 + Lx*0.5 + n;
                if(IntProf > Phase.Fields(i, j, k)[locIndex])
                {
                    for(auto alpha = Phase.Fields(i,j,k).begin();
                             alpha < Phase.Fields(i,j,k).end(); alpha++)
                    {
                        if(alpha->index != locIndex)
                            alpha->value *= 1.0 - IntProf;
                    }
                    Phase.Fields(i, j, k).set(locIndex,IntProf);
                    Phase.Fields(i,j,k).flag = 2;
                }

                i = x0 - Lx*0.5 - n;
                if(IntProf > Phase.Fields(i, j, k)[locIndex])
                {
                    for(auto alpha = Phase.Fields(i,j,k).begin();
                             alpha < Phase.Fields(i,j,k).end(); alpha++)
                    {
                        if(alpha->index != locIndex)
                            alpha->value *= 1.0 - IntProf;
                    }
                    Phase.Fields(i, j, k).set(locIndex,IntProf);
                    Phase.Fields(i,j,k).flag = 2;
                }
            }
        }
    }
    Phase.Finalize(BC);
}

void Initializations::Cylinder(PhaseField& Phase, int PhaseFieldIndex,
                               double Radius, double length, int Axis,
                               double x0, double y0, double z0,
                               Settings& locSettings)
{
    double iWidth = locSettings.iWidth;
    int halfIWidth = iWidth/2;

//    Phase.AddGrainInfo(PhaseIndex);                               //TODO: TO CHECK

    switch(Axis)
    {
        case 0:
        {
    //        SphereFixedIdx(Phase, PhaseIndex, iWidth*0.5, x0-length*0.5+halfIWidth, y0, z0,
    //                locSettings);
    //        SphereFixedIdx(Phase, PhaseIndex, iWidth*0.5, x0+length*0.5-halfIWidth, y0, z0,
    //                locSettings);

            for(int i = x0-length*0.5 + halfIWidth; i < x0+length*0.5 - halfIWidth; i++)
            {
                Disc(Phase, PhaseFieldIndex, Radius, Axis, i,y0,z0, locSettings);
            }
            break;
        }
        case 1:
        {
    //        SphereFixedIdx(Phase, PhaseIndex, iWidth*0.5,
    //                          x0, y0-length*0.5, z0, locSettings);
    //        SphereFixedIdx(Phase, PhaseIndex, iWidth*0.5,
    //                          x0, y0+length*0.5, z0, locSettings);

            for(int j = y0-length*0.5 + halfIWidth; j < y0+length*0.5 - halfIWidth; j++)
            {
                Disc(Phase, PhaseFieldIndex, Radius, Axis, x0,j,z0, locSettings);
            }
            break;
        }
        case 2:
        {
    //        SphereFixedIdx(Phase, PhaseIndex, iWidth*0.5,
    //                          x0, y0, z0-length*0.5, locSettings);
    //        SphereFixedIdx(Phase, PhaseIndex, iWidth*0.5,
    //                          x0, y0, z0+length*0.5, locSettings);

            for(int k = z0-length*0.5 + halfIWidth; k < z0+length*0.5 - halfIWidth; k++)
            {
                Disc(Phase, PhaseFieldIndex, Radius, Axis, x0, y0, k, locSettings);
            }

            break;
        }
        default:
        {
            stringstream message;
            message<<"Axis = "<<Axis<<"! Choose a value between 0 and 2!\n";
            Info::WriteExit(message.str(), "Initialisations", "Cylinder()");
            exit(13);
        }
    }
}

void Initializations::Disc(PhaseField& Phase, int PhaseFieldIndex, double Radius,
         int NormalAxis, double x0, double y0, double z0, Settings& locSettings)
{
    double iWidth = locSettings.iWidth;

//    int locIndex = Phase.AddGrainInfo(PhaseIndex);

    switch(NormalAxis)
    {
        case 0:
        {
            for(int j = y0-Radius-iWidth/2-1; j < y0+Radius+iWidth/2+1; ++j)
            for(int k = z0-Radius-iWidth/2-1; k < z0+Radius+iWidth/2+1; ++k)
            {
                double rad = sqrt((j-y0)*(j-y0)+(k-z0)*(k-z0));
                if (rad < Radius - iWidth*0.5)
                {
                    Phase.Fields(x0, j, k).clear();
                    Phase.Fields(x0, j, k).set(PhaseFieldIndex, 1.0);
                }
                else if (rad < Radius + iWidth*0.5)
                {
                    double IntProf = 0.5 - 0.5*sin(Pi*(rad - Radius)/iWidth);
                    if(IntProf > Phase.Fields(x0, j, k)[PhaseFieldIndex])
                    {
                        for(auto alpha = Phase.Fields(x0,j,k).begin();
                                 alpha < Phase.Fields(x0,j,k).end(); alpha++)
                        {
                            alpha->value *= 1.0 - IntProf;
                        }
                        Phase.Fields(x0, j, k).set(PhaseFieldIndex, IntProf);
                        Phase.Fields(x0, j, k).flag = 2;
                    }
                }
            }
            break;
        }
        case 1:
        {
            for(int i = x0-Radius-iWidth/2-1; i < x0+Radius+iWidth/2+1; ++i)
            for(int k = z0-Radius-iWidth/2-1; k < z0+Radius+iWidth/2+1; ++k)
            {
                double rad = sqrt((i-x0)*(i-x0)+(k-z0)*(k-z0));
                if (rad < Radius - iWidth*0.5)
                {
                    Phase.Fields(i, y0, k).clear();
                    Phase.Fields(i, y0, k).set(PhaseFieldIndex, 1.0);
                }
                else if (rad < Radius + iWidth*0.5)
                {
                    double IntProf = 0.5 - 0.5*sin(Pi*(rad - Radius)/iWidth);
                    if(IntProf > Phase.Fields(i, y0, k)[PhaseFieldIndex])
                    {
                        for(auto alpha = Phase.Fields(i,y0,k).begin();
                                 alpha < Phase.Fields(i,y0,k).end(); alpha++)
                        {
                            alpha->value *= 1.0 - IntProf;
                        }
                        Phase.Fields(i, y0, k).set(PhaseFieldIndex, IntProf);
                        Phase.Fields(i, y0, k).flag = 2;
                    }
                }
            }
            break;
        }
        case 2:
        {
            for(int i = x0-Radius-iWidth/2-1; i < x0+Radius+iWidth/2+1; ++i)
            for(int j = y0-Radius-iWidth/2-1; j < y0+Radius+iWidth/2+1; ++j)
            {
                double rad = sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0));
                if (rad < Radius - iWidth*0.5)
                {
                    Phase.Fields(i, j, z0).clear();
                    Phase.Fields(i, j, z0).set(PhaseFieldIndex, 1.0);
                }
                else if (rad < Radius + iWidth*0.5)
                {
                    double IntProf = 0.5 - 0.5*sin(Pi*(rad - Radius)/iWidth);
                    if(IntProf > Phase.Fields(i, j, z0)[PhaseFieldIndex])
                    {
                        for(auto alpha = Phase.Fields(i,j,z0).begin();
                                 alpha < Phase.Fields(i,j,z0).end(); alpha++)
                        {
                            alpha->value *= 1.0 - IntProf;
                        }
                        Phase.Fields(i, j, z0).set(PhaseFieldIndex, IntProf);
                        Phase.Fields(i, j, z0).flag = 2;
                    }
                }
            }
            break;
        }
        default:
        {
            stringstream message;
            message << "NormalAxis = " << NormalAxis
                    << "! Choose a value between 0 and 2!\n";
            Info::WriteExit(message.str(), "Initializations", "Disc()");
            exit(13);
        }
    }
}

void Initializations::SphereFixedIdx(PhaseField& Phase, int PhaseIndex,
                             double Radius, double x0, double y0, double z0,
                             BoundaryConditions& BC, Settings& locSettings)
{
    double iWidth = locSettings.iWidth;
    int locIndex = PhaseIndex;
    Phase.FieldsStatistics[locIndex].Exist = true;
    Phase.FieldsStatistics[locIndex].Stage = 0;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        double rad = sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0)+(k-z0)*(k-z0));
        if (rad < Radius - iWidth*0.5)
        {
            Phase.Fields(i, j, k).clear();
            Phase.Fields(i, j, k).set(locIndex, 1.0);
        }
        else if (rad < Radius + iWidth*0.5)
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(rad - Radius)/iWidth);
            if(IntProf > Phase.Fields(i, j, z0)[PhaseIndex])
            {
                for(auto alpha = Phase.Fields(i,j,k).begin();
                         alpha < Phase.Fields(i,j,k).end(); alpha++)
                {
                    alpha->value *= 1.0 - IntProf;
                }
                Phase.Fields(i, j, k).set(locIndex, IntProf);
                Phase.Fields(i,j,k).flag = 2;
            }
        }
    }
    STORAGE_LOOP_END

    Phase.Finalize(BC);
}
/*
int Initializations::Wall(PhaseField& Phase, int PhaseIndex, int WallWidth,
                          BoundaryConditions& BC, Settings& locSettings)
{
    const int Nx = locSettings.Nx;
    const int Ny = locSettings.Ny;

    int locIndex = Phase.AddGrainInfo(PhaseIndex);

    for(int i = 0; i < Nx; ++i)
    for(int j = 0; j < Ny; ++j)
    for(int k = 0; k < WallWidth; ++k)
    {
        Phase.Fields(i, j, k).clear();
        Phase.Fields(i, j, k).set(locIndex, 1.0);
    }
    Phase.Finalize(BC);

    return locIndex;
}*/

int Initializations::TwoDimEBSD(std::string filename, std::vector<int> columns,
        PhaseField& Phase, BoundaryConditions& BC, Settings& locSettings)
{
    // Method to read in 2D EBSD data from input file with the following format:
    //    x-coord     y-coordinate     phase-field index (this line is not required in the file)
    // "    0            0                    1       "
    // "    0.1          0                    2       "
    // "                     etc                      "

    // The columns are seperated by whitespaces, the exact column numbers can be given
    // via the columns vector using the following correlations:
    // x coordinate -> columns[0]
    // y coordinate -> columns[1]
    // phase field index -> columns[2]
    //
    // Example call:
    //    Initializations::TwoDimEBSDWithOrientations("EBSDmap.dat",
    //                      {4, 5, 6, 1, 2, 3}, Phi, EP, BC, OPSettings);
    //
    // Note: For hexagonal grids run /scripts/EBSDHexRegParser.m first (Matlab)

    if(columns.size() != 3)
    {
        columns.clear(); columns = {1, 2, 3};
        Info::WriteWarning("Column list incomplete. Set to default {1, 2, 3}.", "Initializations", "TwoDimEBSD");
    }

    if(locSettings.Nz != 1)
    {
        Info::WriteExit("Only 2d supported", "Initializations", "TwoDimEBSD");
        exit(1);
    }

    ifstream inp(filename.c_str(), ios::out);

    if (!inp)
    {
        Info::WriteExit("File " + filename + " could not be opened",
                "Initializations", "TwoDimEBSD()");
        exit(1);
    };

    Info::WriteBlankLine();
    Info::WriteLineInsert("EBSD reader");
    Info::WriteStandard("Source", filename);
    Info::WriteStandard("Read x coordinates from column", std::to_string(columns[0]));
    Info::WriteStandard("Read y coordinates from column", std::to_string(columns[1]));
    Info::WriteStandard("Read phase indices from column", std::to_string(columns[2]));

    // Count number of lines
    int nol = 0;
    std::string line;
    while(getline (inp,line)) nol++;
    inp.clear(); // back to begin
    inp.seekg(0, ios::beg);
    Info::WriteStandard("Number of lines", std::to_string(nol));

    // Read input file

    vector<int> phaseindex;
    vector<int> individualphaseindices;
    vector<double> xcoord;
    vector<double> individualxcoord;
    vector<double> ycoord;
    vector<double> individualycoord;

    while(true)
    {
        std::string line;
        std::getline( inp, line );
        if( !inp ) break;
        std::istringstream iline( line );
//        if (line.size() > 0) cout << line << endl;

        int spos = 1;

        while(true)
        {
            std::string str;
            std::getline(iline, str, ' ' );

            if(spos == columns[0] and str.size() > 0) // Read x coordinates
            {
                double xx = std::stod(str);

                // Identify individual number of xcoordinates
                if(std::find(individualxcoord.begin(), individualxcoord.end(), xx) == individualxcoord.end())
                {
                    individualxcoord.push_back(xx);
    //                cout << xx << endl;
                }
                xcoord.push_back(xx);
            }

            if(spos == columns[1] and str.size() > 0)  // Read y coordinates
            {
                double yy = std::stod(str);
                // Identify individual number of ycoordinates
                if(std::find(individualycoord.begin(), individualycoord.end(), yy) == individualycoord.end())
                {
                    individualycoord.push_back(yy);
                }
                ycoord.push_back(yy);
            }

            if(spos == columns[2] and str.size() > 0) // Read phase field index
            {
                int lindex = std::stoi(str);
                // Identify individual number of phases
                if (std::find(individualphaseindices.begin(), individualphaseindices.end(), lindex) == individualphaseindices.end())
                {
                    individualphaseindices.push_back(lindex);
                }
                phaseindex.push_back(lindex);
            }
            if (str.size() > 0) spos++;
            if( !iline ) break;
        }
    }

    int individualgrainnum = individualphaseindices.size();
    Info::WriteStandard("Number of individual phases", std::to_string(individualgrainnum));
    int maxphaseindex = *max_element(individualphaseindices.begin(), individualphaseindices.end());
    Info::WriteStandard("Largest phase-field index", std::to_string(maxphaseindex));
    int pfsize = phaseindex.size();
    Info::WriteStandard("Phasefield storage size", std::to_string(pfsize));
    int NxEBSD = individualxcoord.size();
    Info::WriteStandard("Number of (individual) x coordinates", std::to_string(NxEBSD));
    double maxX = xcoord.back();
    Info::WriteStandard("Biggest X Coord", std::to_string(maxX));
    double distX = maxX/double(NxEBSD);
    Info::WriteStandard("Rastering X", std::to_string(distX));
    int NyEBSD = individualycoord.size();
    Info::WriteStandard("Number of (individual) y coordinates", std::to_string(NyEBSD));
    double maxY = ycoord.back();
    Info::WriteStandard("Biggest Y Coord", std::to_string(maxY));
    double distY = maxY/double(NyEBSD);
    Info::WriteStandard("Rastering Y", std::to_string(distY));

//    for(size_t i = 0 ; i < phaseindex.size(); ++i) cout << phaseindex[i] << " ";
//    cout << endl;
//
//    for(size_t i = 0 ; i < xcoord.size(); ++i) cout << xcoord[i] << " ";
//    cout << endl;
//
    stringstream outPhaseIndeces;
    for(size_t i = 0 ; i < individualphaseindices.size(); ++i) outPhaseIndeces << individualphaseindices[i] << " ";
    Info::WriteSimple(outPhaseIndeces.str());
    Info::WriteBlankLine();

    const int Nx = locSettings.Nx;
    const int Ny = locSettings.Ny;

    for (int idx = 0; idx < maxphaseindex+1; idx++)
    {
        if(std::find(individualphaseindices.begin(), individualphaseindices.end(), idx) != individualphaseindices.end())
        {
            int locIndex = Phase.AddGrainInfo(0) + idx%locSettings.Nphases; // eigene zv

            Info::WriteStandard("Add grain/phase", std::to_string(locIndex));

            int incx = 1; // Skip points in x-direction
            int incy = 1; // Skip points in y-direction

            int itery = 1;
            int colstep = 0;
            for(int j = 0; j < Ny; ++j)
            {
                int iterx = 1;

                for(int i = 0; i < Nx; ++i)
                {
                    if (idx + 1 == phaseindex[-1 + iterx + colstep])
    //                if ((locIndex)%locSettings.Nphases + 1 + idx == phaseindex[-1 + iterx + colstep])
                    {
    //                    cout << i << " " << j << " " << -1 + iterx + (itery-1)*Nx << endl;
                        Phase.Fields(i, j, 0).clear();
                        Phase.Fields(i, j, 0).flag = 2;
                        Phase.Fields(i, j, 0).set(locIndex, 1.0);
                    }
                    iterx += incx;
                }
                if(incx == 1 and incy == 1)
                {
                    colstep += Nx;
                }
                else
                {
                // Odd columns
                    if(j%2 == 1) colstep += Nx;
                    else colstep += (Nx-1);
                }
                itery += incy;
            }
        }
    }

    Phase.Finalize(BC);
    Info::WriteStandard("EBSD", "done");
    Info::WriteLine();

    return 0;
}

int Initializations::TwoDimEBSDWithOrientations(std::string filename, std::vector<int> columns, std::string anglerepresentation,
        PhaseField& Phase, ElasticProperties& EP, Orientations& OR, BoundaryConditions& BC, Settings& locSettings)
{
    // Method to read in 2D EBSD data from input file with the following format:
    //    x-coord     y-coordinate     phase-field index (this line is not required in the file)
    // "    0            0                    1       "
    // "    0.1          0                    2       "
    // "                     etc                      "

    // The columns are seperated by whitespaces, the exact column numbers can be given
    // via the columns vector using the following correlations:
    // x coordinate -> columns[0]
    // y coordinate -> columns[1]
    // phase field index -> columns[2]
    // Euler angle 1 -> columns[3]
    // Euler angle 2 -> columns[4]
    // Euler angle 3 -> columns[5]
    //
    // Example call:
    //    Initializations::TwoDimEBSDWithOrientations("EBSDmap.dat",
    //                      {4, 5, 6, 1, 2, 3}, Phi, EP, BC, OPSettings);
    //
    // Note: For hexagonal grids run /scripts/EBSDHexRegParser.m first (Matlab)

    const int Nx = locSettings.Nx;
    const int Ny = locSettings.Ny;
    const int Nz = locSettings.Nz;

    if(columns.size() != 6)
    {
        columns.clear(); columns = {1, 2, 3, 4, 5, 6};
        Info::WriteWarning("Column list incomplete. Set to default {1, 2, 3, 4, 5, 6}.",
                "Initializations", "TwoDimEBSD");
    }

    if(locSettings.Nz != 1)
    {
        Info::WriteExit("Only 2d supported", "Initializations", "TwoDimEBSD");
        exit(1);
    }

    ifstream inp(filename.c_str(), ios::out);

    if (!inp)
    {
        Info::WriteExit("File " + filename + " could not be opened",
                "Initializations", "TwoDimEBSD()");
        exit(1);
    };

    Info::WriteBlankLine();
    Info::WriteLineInsert("EBSD reader");
    Info::WriteStandard("Source", filename);
    Info::WriteStandard("Read x coordinates from column", std::to_string(columns[0]));
    Info::WriteStandard("Read y coordinates from column", std::to_string(columns[1]));
    Info::WriteStandard("Read phase indices from column", std::to_string(columns[2]));
    Info::WriteStandard("Read Euler angle 1 from column", std::to_string(columns[3]));
    Info::WriteStandard("Read Euler angle 2 from column", std::to_string(columns[4]));
    Info::WriteStandard("Read Euler angle 3 from column", std::to_string(columns[5]));

    // Count number of lines
    int nol = 0;
    std::string line;
    while(getline (inp,line)) nol++;
    inp.clear(); // back to begin
    inp.seekg(0, ios::beg);
    Info::WriteStandard("Number of lines", std::to_string(nol));

    // Read input file

    vector<int> phaseindex;
    vector<int> individualphaseindices;
    vector<double> xcoord;
    vector<double> individualxcoord;
    vector<double> ycoord;
    vector<double> individualycoord;
    vector<EulerAngles> Eang;

    while(true)
    {
        std::string line;
        std::getline( inp, line );
        if( !inp ) break;
        std::istringstream iline( line );
//        if (line.size() > 0) cout << line << endl;

        EulerAngles tempEang;
        tempEang.set_to_zero();
        int spos = 1;

        while(true)
        {
            std::string str;
            std::getline(iline, str, ' ' );

            if(spos == columns[0] and str.size() > 0) // Read x coordinates
            {
                double xx = std::stod(str);

                // Identify individual number of xcoordinates
                if(std::find(individualxcoord.begin(), individualxcoord.end(), xx) == individualxcoord.end())
                {
                    individualxcoord.push_back(xx);
    //                cout << xx << endl;
                }
                xcoord.push_back(xx);
            }

            if(spos == columns[1] and str.size() > 0)  // Read y coordinates
            {
                double yy = std::stod(str);
                // Identify individual number of ycoordinates
                if(std::find(individualycoord.begin(), individualycoord.end(), yy) == individualycoord.end())
                {
                    individualycoord.push_back(yy);
                }
                ycoord.push_back(yy);
            }

            if(spos == columns[2] and str.size() > 0) // Read phase field index
            {
                int lindex = std::stoi(str);
                // Identify individual number of phases
                if (std::find(individualphaseindices.begin(), individualphaseindices.end(), lindex) == individualphaseindices.end())
                {
                    individualphaseindices.push_back(lindex);
                }
                phaseindex.push_back(lindex);
            }

            double anglefactor = 180.0/Pi;
            if(!anglerepresentation.compare("degree") or !anglerepresentation.compare("deg")
                    or !anglerepresentation.compare("Degree") or !anglerepresentation.compare("Deg"))
            {
                anglefactor = 1.0;
            }

            if(spos == columns[3] and str.size() > 0) // Read Euler angle 1
            {
                tempEang.Q[0] = std::stod(str)*anglefactor;
            }

            if(spos == columns[4] and str.size() > 0) // Read Euler angle 2
            {
                tempEang.Q[1] = std::stod(str)*anglefactor;
            }

            if(spos == columns[5] and str.size() > 0) // Read Euler angle 3
            {
                tempEang.Q[2] = std::stod(str)*anglefactor;
            }

            if (str.size() > 0) spos++;
            if( !iline ) break;
        }

        tempEang.set_convention(ZXZ);

        Info::WriteWarning("Using ZXZ angle convention", "Initializations", "TwoDimEBSDWithOrientations()");

        tempEang.setTrigonometricFunctions();
        Eang.push_back(tempEang);
    }

    int individualgrainnum = individualphaseindices.size();
    Info::WriteStandard("Number of individual phases", std::to_string(individualgrainnum));
    int pfsize = phaseindex.size();
    Info::WriteStandard("Phasefield storage size", std::to_string(pfsize));
    int NxEBSD = individualxcoord.size();
    Info::WriteStandard("Number of (individual) x coordinates", std::to_string(NxEBSD));
    double maxX = xcoord.back();
    Info::WriteStandard("Biggest X Coord", std::to_string(maxX));
    double distX = maxX/double(NxEBSD);
    Info::WriteStandard("Rastering X", std::to_string(distX));
    int NyEBSD = individualycoord.size();
    Info::WriteStandard("Number of (individual) y coordinates", std::to_string(NyEBSD));
    double maxY = ycoord.back();
    Info::WriteStandard("Biggest Y Coord", std::to_string(maxY));
    double distY = maxY/double(NyEBSD);
    Info::WriteStandard("Rastering Y", std::to_string(distY));

    Storage3D <EulerAngles, 0> StorageEulerAngles;
    StorageEulerAngles.Allocate(Nx, Ny, 1, 0);

    for (int idx = 0; idx < individualgrainnum; idx++)
    {
        int locIndex = Phase.AddGrainInfo(0) + idx%locSettings.Nphases;

        Info::WriteStandard("Add grain/phase", std::to_string(locIndex));

        int incx = 1; // Skip points in x-direction
        int incy = 1; // Skip points in y-direction

        int itery = 1;
        int colstep = 0;
        for(int j = 0; j < Ny; ++j)
        {
            int iterx = 1;

            for(int i = 0; i < Nx; ++i)
            {
                // SetInitialOrientations
                EulerAngles EangLocal = Eang[-1 + iterx + colstep];
                OR.Quaternions(i,j,0) = EangLocal.getQuaternion();
                StorageEulerAngles(i,j,0) = EangLocal;

                // Set phase fields
                if (idx + 1 == phaseindex[-1 + iterx + colstep])
//                if ((locIndex)%locSettings.Nphases + 1 + idx == phaseindex[-1 + iterx + colstep])
                {
//                    cout << i << " " << j << " " << -1 + iterx + (itery-1)*Nx << endl;
                    Phase.Fields(i, j, 0).clear();
                    Phase.Fields(i, j, 0).flag = 2;
                    Phase.Fields(i, j, 0).set(locIndex, 1.0);
                }
                iterx += incx;
            }
            if(incx == 1 and incy == 1)
            {
                colstep += Nx;
            }
            else
            {
            // Odd columns
                if(j%2 == 1) colstep += Nx;
                else colstep += (Nx-1);
            }
            itery += incy;
        }
    }

    // Create output to VTK

    stringstream outbuffer;

    outbuffer << "# vtk DataFile Version 3.0\n";
    outbuffer << "InitialEulerAngles\n";
    outbuffer << "ASCII\n";
    outbuffer << "DATASET STRUCTURED_GRID\n";
    outbuffer << "DIMENSIONS " << Nx << " " << Ny << " " << 1 << "\n";
    outbuffer << "POINTS " <<  Nx*Ny << " int\n";

    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        outbuffer << i << " " << j << " " << k << "\n";
    }
    outbuffer << " \n";
    outbuffer << "POINT_DATA " << Nx*Ny << " \n";

    outbuffer << "SCALARS Eang_" << 1 << " double\n";
    outbuffer << "LOOKUP_TABLE default\n";
    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        outbuffer << StorageEulerAngles(i,j,k).Q[0] << " ";
    }
    outbuffer << " \n";
    outbuffer << "SCALARS Eang_" << 2 << " double\n";
    outbuffer << "LOOKUP_TABLE default\n";
    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        outbuffer << StorageEulerAngles(i,j,k).Q[1] << " ";
    }
    outbuffer << " \n";
    outbuffer << "SCALARS Eang_" << 3 << " double\n";
    outbuffer << "LOOKUP_TABLE default\n";
    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        outbuffer << StorageEulerAngles(i,j,k).Q[2] << " ";
    }

    string FileName = VTKDir + "InitializedEulerAngles.vtk";

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbuffer.rdbuf();
    vtk_file.close();

    // Write table of orientations

    outbuffer.str("");

    for(size_t i = 0 ; i < Eang.size(); ++i)
    {
        outbuffer << i << "  " << Eang[i].Q[0] << "  "
                              << Eang[i].Q[1] << "  "
                              << Eang[i].Q[2] << endl;
    }
    FileName = "EBSDorientations.dat";

    ofstream orientations_file(FileName.c_str());
    orientations_file << outbuffer.rdbuf();
    orientations_file.close();

    // End write Euler angles

    Phase.Finalize(BC);
    Info::WriteStandard("EBSD", "done");
    Info::WriteLine();

    return 0;
}

int Initializations::SphereInGrain(PhaseField& Phase, int ParrentPhaseFieldIndex,
                 int PhaseIndex, double Radius, double x0, double y0, double z0,
                                  BoundaryConditions& BC, Settings& locSettings)
{
    double iWidth = locSettings.iWidth;
    int locIndex = Phase.AddGrainInfo(PhaseIndex);
    cout<<locIndex<<endl;
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        double ParrentPhaseAmount = Phase.Fields(i, j, k)[ParrentPhaseFieldIndex];

        double rad = sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0)+(k-z0)*(k-z0));
        if (rad < (Radius - iWidth*0.5))
        {
            Phase.Fields(i,j,k).set(ParrentPhaseFieldIndex, 0.0);
            Phase.Fields(i,j,k).set(locIndex, ParrentPhaseAmount);
            Phase.Fields(i,j,k).flag = 2;
        }
        else if (rad < (Radius + iWidth*0.5))
        {
            double Profile = (0.5 - 0.5*sin(Pi*(rad - Radius)/iWidth))*ParrentPhaseAmount;

            Phase.Fields(i, j, k).set(ParrentPhaseFieldIndex, ParrentPhaseAmount-Profile);
            Phase.Fields(i, j, k).set(locIndex, Profile);
            Phase.Fields(i,j,k).flag = 2;
        }
    }
    STORAGE_LOOP_END

    Phase.Finalize(BC);

    return locIndex;
}

int  Initializations::FillGrainWithRandomSpheres(PhaseField& Phase,
                              int ParrentPhaseFieldIndex, int SpheresPhaseIndex,
               double x0, double y0, double z0, double Lx, double Ly, double Lz,
        double MinR, double MaxR, BoundaryConditions& BC, Settings& locSettings)
{
    vector<dVector3> spheres;
    vector<double> Radii;

    int counter = 0;
    int globalCounter = 0;
    const int MaxNumberOfTries = 10*(Lx-x0)*(Ly-y0)*(Lz-z0);
//    const int MaxNumberOfTries = 10000000;

    default_random_engine generator;
    uniform_int_distribution<int> distributionX(x0, x0 + Lx-1);
    uniform_int_distribution<int> distributionY(y0, y0 + Ly-1);
    uniform_int_distribution<int> distributionZ(z0, z0 + Lz-1);
    uniform_real_distribution<double> distributionR(MinR, MaxR);
    while(counter < MaxNumberOfTries)
    {
        dVector3 coordinates;
        coordinates[0] = distributionX(generator);
        coordinates[1] = distributionY(generator);
        coordinates[2] = distributionZ(generator);
//        cout<<x<<","<<y<<","<<z<<endl;

        bool overlapping = false;
        double Radius = distributionR(generator);

        for(size_t n = 0; n < spheres.size(); n++)
        {
            double distance = (spheres[n]-coordinates).abs();
            if(distance < Radius+Radii[n] + 10) overlapping = true;
        }

        if(overlapping or
           Phase.Fields(coordinates[0], coordinates[1], coordinates[2])[ParrentPhaseFieldIndex] < 1.0)
        {
            counter++;
        }
        else
        {
            int sphereIndex =
            SphereInGrain(Phase, ParrentPhaseFieldIndex, SpheresPhaseIndex,
                         Radius, coordinates[0], coordinates[1], coordinates[2],
                                                               BC, locSettings);
            spheres.push_back(coordinates);
            Radii.push_back(Radius);
            stringstream message;
            message << "*********************************************************\n";
            message << "after " << counter << " tries sphere nr: "
                    << spheres.size() << " at " << coordinates.print()
                    << " with R = " << Radius<< " and phase field index = "
                    << sphereIndex << " set." << ";\n";
            Phase.PrintPointStatistics(coordinates[0], coordinates[1], coordinates[2]);
            message << "*********************************************************\n";
            Info::WriteSimple(message.str());
            globalCounter += counter;
            counter = 0;
        }
    }
    cout << "after " << globalCounter << " tries " << spheres.size()
         << " spheres could be initialized!\n";

    Phase.Finalize(BC);

    return spheres.size();
}

void RemoveParrentGrain(PhaseField& Phase, int ParrentPhaseFieldIndex,
                        BoundaryConditions& BC, Settings& locSettings)
{
    /*
     * Removes phase number ParrentPhaseFieldIndex and increases phase number 0
     * by the removed amount.
     */
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        double ParrentPhaseAmount = Phase.Fields(i, j, k)[ParrentPhaseFieldIndex];

        if(ParrentPhaseAmount != 0.0)
        {
            Phase.Fields(i, j, k).add(0, ParrentPhaseAmount);
            Phase.Fields(i, j, k).set(ParrentPhaseFieldIndex, 0);
        }
    }
    STORAGE_LOOP_END
    Phase.Finalize(BC);
}

int  Initializations::FillGrainWithSpheres(PhaseField& Phase,
                              int ParrentPhaseFieldIndex, int SpheresPhaseIndex,
               double x0, double y0, double z0, double Lx, double Ly, double Lz,
             double SphereRadius, BoundaryConditions& BC, Settings& locSettings)
{
    int SpheresNumber = 0;

    vector<dVector3> spheres;

    int counter = 0;
    int globalCounter = 0;
    const int MaxNumberOfTries = 10*(Lx-x0)*(Ly-y0)*(Lz-z0);
//    const int MaxNumberOfTries = 10000000;

    default_random_engine generator;
    uniform_int_distribution<int> distributionX(x0, x0 + Lx-1);
    uniform_int_distribution<int> distributionY(y0, y0 + Ly-1);
    uniform_int_distribution<int> distributionZ(z0, z0 + Lz-1);

//    random_device randomNumber;

    while(counter < MaxNumberOfTries)
    {
        dVector3 coordinates;
        coordinates[0] = distributionX(generator);
        coordinates[1] = distributionY(generator);
        coordinates[2] = distributionZ(generator);
//        cout<<x<<","<<y<<","<<z<<endl;

        double minDistance = 1000.0;

        for(unsigned int n = 0; n < spheres.size(); n++)
        {
            double distance = (spheres[n]-coordinates).abs();
            if(distance<minDistance) minDistance = distance;
        }

        if(minDistance < 2*SphereRadius or
           Phase.Fields(coordinates[0], coordinates[1], coordinates[2])[ParrentPhaseFieldIndex] < 1.0)
        {
            counter++;
        }
        else
        {
            stringstream message;
            message<<"*********************************************************\n";
            message<<"after "<<counter<<" tries sphere nr: "<<spheres.size()<<" set." <<";\n";
            Phase.PrintPointStatistics(coordinates[0], coordinates[1], coordinates[2]);
            message<<"*********************************************************\n";
            Info::WriteSimple(message.str());
            SphereInGrain(Phase, ParrentPhaseFieldIndex, SpheresPhaseIndex,
                   SphereRadius, coordinates[0], coordinates[1], coordinates[2],
                                                               BC, locSettings);
//            Sphere(Phase, SpheresPhaseIndex,
//                   SphereRadius, coordinates[0], coordinates[1], coordinates[2],
//                                                               BC, locSettings);
            spheres.push_back(coordinates);
            globalCounter += counter;
            counter = 0;
        }
    }

    stringstream message;
    message << "after " << globalCounter << " tries " << spheres.size()
         << " spheres could be initialized!\n";
    Info::WriteSimple(message.str());

    RemoveParrentGrain(Phase, ParrentPhaseFieldIndex, BC, locSettings);

    Phase.Finalize(BC);

    return SpheresNumber;
}

int  Initializations::FillGrainWithSpheres(PhaseField& Phase,
     int ParrentPhaseFieldIndex, int SpheresPhaseIndex, double SphereRadius,
                             BoundaryConditions& BC, Settings& locSettings)
{
    return FillGrainWithSpheres(Phase, ParrentPhaseFieldIndex,
              SpheresPhaseIndex, locSettings.Nx, locSettings.Ny, locSettings.Nz,
                                        1, 1, 1, SphereRadius, BC, locSettings);
}

int Initializations::DiffusionCouple(PhaseField& Phase,
                                 int MajorityPhaseIndex, int MinorityPhaseIndex,
                                 double MinorityPhaseLayerThickness,
                                 BoundaryConditions& BC, Settings& locSettings)
{
    /** This will set up a diffusion couple, where the minority phase occupies
    a layer in the middle of the box (in x direction) with a given thickness.
    The majority phase occupies the space left and right to this layer. That
    way, periodic boundary conditions can be used in all directions.*/

    int Nx = locSettings.Nx;
    double iWidth = locSettings.iWidth;

    int index1 = Phase.AddGrainInfo(MajorityPhaseIndex);
    int index2 = Phase.AddGrainInfo(MinorityPhaseIndex);

    int offset = MinorityPhaseLayerThickness;

    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();
        Phase.Fields(i,j,k).flag = 0;

        if (i < Nx/2 - offset/2 - iWidth/2)
        {
            Phase.Fields(i, j, k).set(index1, 1);   // left phase
        }
        else if (i > Nx/2 - offset/2 + iWidth/2 and i < Nx/2 + offset/2 - iWidth/2)
        {
            Phase.Fields(i, j, k).set(index2, 1);   // middle phase
        }
        else if (i > Nx/2 + offset/2 + iWidth/2)
        {
            Phase.Fields(i, j, k).set(index1, 1);   // right phase
        }
        else if (i >= Nx/2 - offset/2 - iWidth/2 and i <= Nx/2 - offset/2 + iWidth/2)
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(i - (Nx/2 - offset/2))/iWidth);  // left interface
            Phase.Fields(i, j, k).set(index2, 1.0 - IntProf);
            Phase.Fields(i, j, k).set(index1, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
        else if (i >= Nx/2 + offset/2 - iWidth/2 and i <= Nx/2 + offset/2 + iWidth/2)
        {
            double IntProf = 0.5 - 0.5*sin(Pi*(i - (Nx/2 + offset/2))/iWidth);  // right interface
            Phase.Fields(i, j, k).set(index1, 1.0 - IntProf);
            Phase.Fields(i, j, k).set(index2, IntProf);
            Phase.Fields(i,j,k).flag = 2;
        }
        else
        {
            cout << "Initialization of <Diffusion Couple> failed. Exiting now." << endl;
            exit(1);
        }
    }
    STORAGE_LOOP_END

    Phase.Finalize(BC);

    return index1*10+index2;
}
/*
int Initializations::MechanicalWorkBench(PhaseField& Phase, int nrVoidCells, int voidPhIDX, int matrixPhIDX, BoundaryConditions& BC, Settings& locSettings)
{
    //Initialize soft material
    int locIndex = Phase.AddGrainInfo(voidPhIDX);
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,Phase.Fields.Bcells())
    {
        Phase.Fields(i, j, k).clear();
        Phase.Fields(i, j, k).set(locIndex, 1.0);
    }
    STORAGE_LOOP_END
    //Initialize matrix material
    int Nx = locSettings.Nx;
    int Ny = locSettings.Ny;
    int Nz = locSettings.Nz;

    locIndex = Phase.AddGrainInfo(matrixPhIDX);
    OP_STORAGE_LOOP_RANGE(i,j,k,0,nrVoidCells-1,0,Nx,Ny-nrVoidCells,Nz)
    {
        Phase.Fields(i, j, k).clear();
        Phase.Fields(i, j, k).set(locIndex, 1.0);
    }

    Phase.Finalize(BC);
    return locIndex;
}
*/
void Initializations::VoronoiTesselation(PhaseField& Phase, BoundaryConditions& BC,
	Settings& OPSettings, int Ngrains, int MatrixPhase)
{
	// 3D version of Voronoi initialization.
#ifdef NoVoro
	Single(Phase, MatrixPhase, BC, OPSettings);
#else
	Info::WriteLineInsert("Voronoi creator");
	Info::WriteStandard("Number of grains", std::to_string(Ngrains));
	Info::WriteStandard("Thermodynamic phase", std::to_string(MatrixPhase));

	Phase.Clear();

	double rX0 = 0.0;
	double rXend = (OPSettings.Nx - 1);
	double rY0 = 0.0;
	double rYend = (OPSettings.Ny - 1);
	double rZ0 = 0.0;
	double rZend = (OPSettings.Nz - 1);

	int Nx = OPSettings.Nx;
	int Ny = OPSettings.Ny;
	int Nz = OPSettings.Nz;

	bool Xperiodic = false;
	bool Yperiodic = false;
	bool Zperiodic = false;

	if (BC.BC0X == Periodic) Xperiodic = true;
	if (BC.BC0Y == Periodic) Yperiodic = true;
	if (BC.BC0Z == Periodic) Zperiodic = true;

	if (!Xperiodic or !Yperiodic or !Zperiodic)
	{
		Info::WriteExit("Voronoi Tesselation needs periodic boundary conditions",
			"Initializations", "VoronoiTesselation");
		exit(1);
	}

	vector<int> individualphaseindices;
	// Create voro++ container with given boundary conditions
	container con(rX0, rXend, rY0, rYend, rZ0, rZend, Nx, Ny, Nz, Xperiodic, Yperiodic, Zperiodic, 8);

	// Randomly add particles into the voro++ container
	for (int i = 0; i < Ngrains; i++)
	{
		con.put(i, rXend*rnd(), rYend*rnd(), rZend*rnd());
	}

	double rx, ry, rz;
	// Populate phasefield storage with grains (sharp interface)
	for (double z = rZ0; z <= rZend; z++)
		for (double y = rY0; y <= rYend; y++)
			for (double x = rX0; x <= rXend; x++)
			{
				int indexreturn;
				if (con.find_voronoi_cell(x, y, z, rx, ry, rz, indexreturn))
				{
					if (indexreturn >= (int)individualphaseindices.size())
					{
						Phase.AddGrainInfo(MatrixPhase);
						individualphaseindices.push_back(indexreturn);
					}
					Phase.Fields(x, y, z).set(indexreturn, 1.0);
				}
			}
	// Spread interfaces by overlaping the grains near the interface
	STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields, 0)
		if (Phase.Fields(i, j, k).size() == 1)
		{
			int alpha_index = Phase.Fields(i, j, k).begin()->index;

			for (int di = -1; di <= 1; di++)
				for (int dj = -1; dj <= 1; dj++)
					for (int dk = -1; dk <= 1; dk++)
						if (i + di >= 0 and i + di < Nx and
							j + dj >= 0 and j + dj < Ny and
							k + dk >= 0 and k + dk < Nz and
							(di != 0 or dj != 0 or dk != 0))
						{
							bool grain_not_present = true;
							for (auto beta = Phase.Fields(i + di, j + dj, k + dk).begin();
								beta < Phase.Fields(i + di, j + dj, k + dk).end(); ++beta)
								if (beta->index == alpha_index)
								{
									grain_not_present = false;
									break;
								}
							if (grain_not_present)
							{
								for (auto beta = Phase.Fields(i + di, j + dj, k + dk).begin();
									beta < Phase.Fields(i + di, j + dj, k + dk).end(); ++beta)
								{
									Phase.Fields(i, j, k).add(beta->index, 1.0);
								}
								Phase.Fields(i + di, j + dj, k + dk).add(alpha_index, 1.0);
							}
						}
		}
	STORAGE_LOOP_END

		Phase.Finalize(BC);

	mt19937_64 OrientGenerator1(45);
	mt19937_64 OrientGenerator2(697);
	mt19937_64 OrientGenerator3(255);

	uniform_int_distribution <int> Q1Distribution(0, 360);
	uniform_int_distribution <int> Q2Distribution(0, 180);
	uniform_int_distribution <int> Q3Distribution(0, 360);

	for (unsigned int n = 0; n < Phase.FieldsStatistics.size(); n++)
	{
		double Q1 = Q1Distribution(OrientGenerator1) * Pi / 180.0;
		double Q2 = Q2Distribution(OrientGenerator2) * Pi / 180.0;
		double Q3 = Q3Distribution(OrientGenerator3) * Pi / 180.0;
		EulerAngles locAngles({ Q1, Q2, Q3 }, XYZ);
		Phase.FieldsStatistics[n].Orientation = locAngles.getQuaternion();
	}

	Info::WriteLine();
#endif
}

void Initializations::ReadSubset(PhaseField& Phase, string PFFileName, string GSFileName,
                                             int offsetX,
                                             int offsetY,
                                             int offsetZ,
                                             const BoundaryConditions& BC,
                                             int inpNx,
                                             int inpNy,
                                             int inpNz)
{
    ReadSubsetGeneral(Phase, PFFileName, GSFileName, offsetX, offsetY, offsetZ,
            inpNx-offsetX, inpNy-offsetY, inpNz-offsetZ, inpNx, inpNy, inpNz,
            0, 0, 0, {false}, BC);
}

void Initializations::ReadSubsetGeneral(PhaseField& Phase, string PFFileName, string GSFileName,
                                             int offsetInpX,
                                             int offsetInpY,
                                             int offsetInpZ,
                                             int sizeInpX,
                                             int sizeInpY,
                                             int sizeInpZ,
                                             int totalSizeInpNx,
                                             int totalSizeInpNy,
                                             int totalSizeInpNz,
                                             int offsetLocalX,
                                             int offsetLocalY,
                                             int offsetLocalZ,
                                             std::initializer_list<bool> newGrain,
                                             const BoundaryConditions& BC)
{
    fstream inp(PFFileName.c_str(), ios::in | ios::binary);

    GrainInfo FieldsStatisticsInp;
    FieldsStatisticsInp.Read(GSFileName);
    int FSInpSize = FieldsStatisticsInp.size();
    int allocateNewGrains = Phase.FieldsStatistics.size();
    vector<int> targetPF(FSInpSize);                                            //Target phase field for every phase field in the source file

    for (int i = 0; i < FSInpSize; i++)
    {
        if ((unsigned)i > newGrain.size() -1)
        {
            Phase.FieldsStatistics.add_grain(FieldsStatisticsInp[i].Phase);
            targetPF[i] = Phase.FieldsStatistics.size()-1;                            //add the PF index of the current Grain to target PF
            allocateNewGrains++;
        }
        else if (newGrain.begin()[i])
        {
            Phase.FieldsStatistics.add_grain(FieldsStatisticsInp[i].Phase);
            targetPF[i] = Phase.FieldsStatistics.size()-1;                            //add the PF index of the current Grain to target PF
            allocateNewGrains++;
        }
        else
        {
            if (!Phase.FieldsStatistics[i].Exist or 
                 Phase.FieldsStatistics[i].Phase != FieldsStatisticsInp[i].Phase)
            {
                stringstream message;
                message << "Target phase-field does not exist!\n"
                        << "Grain No.: " << i << "\n";
                Info::WriteExit(message.str(),
                        "Initializations", "ReadSubsetGeneral()");
                exit(1);
            }
            targetPF[i] = i;                                                    //do not change the PF index
        }
    }

    if (!inp)
    {
        Info::WriteExit(PFFileName + " could not be opened",
                "Initializations", "ReadSubsetGeneral()");
        exit(1);
    };

    int locNx = Phase.Nx;
    int locNy = Phase.Ny;
    int locNz = Phase.Nz;

    inp.read(reinterpret_cast<char*>(&locNx), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNy), sizeof(int));
    inp.read(reinterpret_cast<char*>(&locNz), sizeof(int));
  
    if(Phase.Nx - offsetLocalX < sizeInpX or
       Phase.Ny - offsetLocalY < sizeInpY or
       Phase.Nz - offsetLocalZ < sizeInpZ)
    {
        stringstream message;
        message << "Inconsistent system dimensions!\n"
                << "Input data dimensions: (" << locNx << ", "
                << locNy << ", " << locNz << ") grid points.\n"
                << "Local data dimensions: (" 
                << Phase.Nx << ", " << Phase.Ny << ", " << Phase.Nz 
                << ") grid points.\n"
                << "Local data offset relative to the input data: ("
                << offsetInpX << ", " << offsetInpY << ", " << offsetInpZ
                << ") grid points.\n";
        Info::WriteExit(message.str(), "Initializations", "ReadSubsetGeneral()");
        exit(1);
    }

    for(int i = 0; i < locNx; i++)
    for(int j = 0; j < locNy; j++)
    for(int k = 0; k < locNz; k++)
    {
        Node locPF;
        int   num = 0;
        int   idx = 0;
        double val = 0.0;

        inp.read(reinterpret_cast<char*>(&num), sizeof(int));

        for(int n = 0; n < num; n++)
        {
            inp.read(reinterpret_cast<char*>(&idx), sizeof(int));
            inp.read(reinterpret_cast<char*>(&val), sizeof(double));
            locPF.set(targetPF[idx], val);
        }

        if((i >= offsetInpX and j >= offsetInpY and k >= offsetInpZ) and
           (i <= offsetInpX + sizeInpX and 
            j <= offsetInpY + sizeInpY and
            k <= offsetInpZ + sizeInpZ))
        {
            int x = i - offsetInpX + offsetLocalX;
            int y = j - offsetInpY + offsetLocalY;
            int z = k - offsetInpZ + offsetLocalZ;

            Phase.Fields(x,y,z) = locPF;
            Phase.Fields(x,y,z).flag = Phase.Fields(x,y,z).finalize();
        }
    }

    inp.close();
    Phase.Finalize(BC);
    Info::WriteStandard("Initializations", "Binary input loaded");
}
}