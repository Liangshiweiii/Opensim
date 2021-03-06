#include "Settings.h"
#include "SystemDimensions.h"
#include "Info.h"
#include "Tools/UserInterface.h"
#include "Tools/OPObjects.h"
#include "Tools/Includes.h"

namespace opensim
{

using namespace std;

void SystemDimensions::Initialize(void)
{
    thisclassname = "SystemDimensions";

    //DefaultInputFileName = ProjectInputDir + "ProjectInput.opi";
    
    Nx = 1;
    Ny = 1;
    Nz = 1;

    dNx = 1;
    dNy = 1;
    dNz = 1;
    
    newNx = Nx;
    newNy = Ny;
    newNz = Nz;
    
    OffsetX = 0;
    OffsetY = 0;
    OffsetZ = 0;
    
    maxNx = Nx;
    maxNy = Ny;
    maxNz = Nz;
    
    dx = 1;
    iWidth = 1;

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void SystemDimensions::ReadInput()
{
    ReadInput(DefaultInputFileName);
}

void SystemDimensions::ReadInput(string InputFileName)
{    
    if(thisclassname == "")
    {
        Info::WriteExit("Attempt to use uninitialized object", "SystemDimensions", "ReadInput()");
        exit(1);
    }
    
    fstream inp(InputFileName.c_str(), ios::in);

    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput");
        exit(1);
    };
    
    Info::WriteLine();
    Info::WriteLineInsert("System Dimensions");
    Info::WriteStandard("Source", InputFileName);

    //int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
	int moduleLocation = UserInterface::FindModuleLocation(inp, "Settings");

    Nx = UserInterface::ReadParameterI(inp, moduleLocation, string("Nx"));
    Ny = UserInterface::ReadParameterI(inp, moduleLocation, string("Ny"));
    Nz = UserInterface::ReadParameterI(inp, moduleLocation, string("Nz"));
    
    if(Nx == 0) // reduced X dimension
	{
    	Nx  = 1;
    	dNx = 0;
	}
    if(Ny == 0) // reduced Y dimension
	{
		Ny  = 1;
		dNy = 0;
	}
    if(Nz == 0)  // reduced Z dimension
	{
		Nz  = 1;
		dNz = 0;
	}
    newNx = Nx;
    newNy = Ny;
    newNz = Nz;
    
    maxNx = Nx;
    maxNy = Ny;
    maxNz = Nz;
    
    History.push_back(std::array<int, 4>({0, Nx, Ny, Nz}));
    
    dx = UserInterface::ReadParameterD(inp, moduleLocation, string("dx"));
    iWidth = UserInterface::ReadParameterD(inp, moduleLocation, string("IWidth"));

    inp.close();
    Info::WriteLine();
}

void SystemDimensions::Resize(int size_x, int size_y, int size_z, int tStep)
{
    newNx = size_x;
    newNy = size_y;
    newNz = size_z;
    
    if(RemeshingAllowed and (Nx != newNx or Ny != newNy or Nz != newNz))
    {
        RemeshAll();
        
        Nx = newNx;
        Ny = newNy;
        Nz = newNz;
        
        maxNx = max(Nx, maxNx);
        maxNy = max(Ny, maxNy);
        maxNz = max(Nz, maxNz);
        
        History.push_back(std::array<int, 4>({tStep, Nx, Ny, Nz}));
    }
}

void SystemDimensions::Read(int tStep)
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"SystemDimensions_", tStep, ".dat");

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be opened", thisclassname, "Read");
        exit(1);
    };
    string tmp;
    getline(inp, tmp); // discards the headline
    
    while(!inp.eof())
    {
        int loctStep;
        
        inp >> loctStep >> Nx >> Ny >> Nz;
        
        History.push_back(std::array<int, 4>({loctStep, Nx, Ny, Nz}));
        
        maxNx = max(Nx, maxNx);
        maxNy = max(Ny, maxNy);
        maxNz = max(Nz, maxNz);
    }
}

void SystemDimensions::Write(int tStep)
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"SystemDimensions_", tStep, ".dat");

    ofstream out(FileName.c_str(), ios::out);

    if (!out)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be created", thisclassname, "Write");
        exit(1);
    };
    
    out << "tStep" << "\t" << "Nx" << "\t" << "Ny" << "\t" << "Nz" << endl;
    
    for(unsigned int n = 0; n < History.size(); n++)
    {
        out << History[n][0] << "\t" << History[n][1] << "\t" << History[n][2] << "\t" << History[n][3] << endl;
    }    
}

void SystemDimensions::AddForRemesh(OPObject& Obj)
{
    ObjectsToRemesh.push_back(&Obj);
    NamesOfObjectsToRemesh.push_back(Obj.thisclassname);
}

void SystemDimensions::RemeshAll(void)
{
    for(unsigned int n = 0; n < ObjectsToRemesh.size(); n++)
    {
        ObjectsToRemesh[n]->Remesh(newNx, newNy, newNz);    
    }
}

void SystemDimensions::Remesh(std::string ObjNameBase, int size_x, int size_y, int size_z)
{
    for(unsigned int n = 0; n < ObjectsToRemesh.size(); n++)
    if(NamesOfObjectsToRemesh[n] == ObjNameBase)
    {
        ObjectsToRemesh[n]->Remesh(size_x, size_y, size_z);    
    }
}    

SystemDimensions& SystemDimensions::operator= (const SystemDimensions& rhs)
{
    if (this != &rhs) // protect against self-assignment
    {
        thisclassname = rhs.thisclassname;
        //DefaultInputFileName = rhs.DefaultInputFileName;
        LengthUnits = rhs.LengthUnits;
        
        Nx = rhs.Nx;
        Ny = rhs.Ny;
        Nz = rhs.Nz;
        
        dNx = rhs.dNx;
        dNy = rhs.dNy;
        dNz = rhs.dNz;

        newNx = rhs.newNx;
        newNy = rhs.newNy;
        newNz = rhs.newNz;
        
        TotalX = rhs.TotalX;
        TotalY = rhs.TotalY;
        TotalZ = rhs.TotalZ;
        
        maxNx = rhs.maxNx;
        maxNy = rhs.maxNy;
        maxNz = rhs.maxNz;        
        
        OffsetX = rhs.OffsetX;
        OffsetY = rhs.OffsetY;
        OffsetZ = rhs.OffsetZ;
        
        dx = rhs.dx;
        iWidth = rhs.iWidth;
        RemeshingAllowed = rhs.RemeshingAllowed;
        
        History = rhs.History;
        ObjectsToRemesh = rhs.ObjectsToRemesh;
        NamesOfObjectsToRemesh = rhs.NamesOfObjectsToRemesh;
    }
    return *this;
}

} //namespace opensim