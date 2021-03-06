#include "Settings.h"
#include "Info.h"
#include "Tools/UserInterface.h"

namespace opensim
{

using namespace std;

/**Directory Structure Settings:*/

void Settings::Initialize(void)
{
    thisclassname = "Settings";

    //DefaultInputFileName = ProjectInputDir + "ProjectInput.opi";

    Info::WriteStartScreen();
    Info::WriteStandard(thisclassname, "Initialized");
}
void Settings::ReadInput()
{
    ReadInput(DefaultInputFileName);
}

void Settings::ReadInput(string InputFileName)
{
    if(thisclassname == "")
    {
        Info::WriteExit("Attempt to use uninitialized object", "Settings", "ReadInput()");
        exit(1);
    }

    fstream inp(InputFileName.c_str(), ios::in);

    if (!inp)
    {
        Info::WriteExit("File " + InputFileName + " could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };

    Info::WriteLine();
    Info::WriteLine();
    Info::WriteStandard("Global project parameters source", InputFileName);
    
	int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
	
    SimulationTitle          = UserInterface::ReadParameterS(inp, moduleLocation, string("SimTtl"));
    Info::WriteLineInsert("Model parameters");
    UnitsOfLength            = UserInterface::ReadParameterF(inp, moduleLocation, string("LUnits"), false, "m");
    UnitsOfMass              = UserInterface::ReadParameterF(inp, moduleLocation, string("MUnits"), false, "kg");
    UnitsOfTime              = UserInterface::ReadParameterF(inp, moduleLocation, string("TUnits"), false, "s");
    UnitsOfEnergy            = UserInterface::ReadParameterF(inp, moduleLocation, string("EUnits"), false, "J");
    nOMP                     = UserInterface::ReadParameterI(inp, moduleLocation, string("nOMP"));

    omp_set_num_threads(nOMP);

    nSteps                   = UserInterface::ReadParameterI(inp, moduleLocation, string("nSteps"));
    dt                       = UserInterface::ReadParameterD(inp, moduleLocation, string("dt"));
    Nx                       = UserInterface::ReadParameterI(inp, moduleLocation, string("Nx"));
    Ny                       = UserInterface::ReadParameterI(inp, moduleLocation, string("Ny"));
    Nz                       = UserInterface::ReadParameterI(inp, moduleLocation, string("Nz"));
    dx                       = UserInterface::ReadParameterD(inp, moduleLocation, string("dx"));
    Info::WriteLineInsert("Phase field parameters");
    iWidth                   = UserInterface::ReadParameterD(inp, moduleLocation, string("IWidth"));
    Info::WriteLineInsert("Additional parameters");
    //T                        = UserInterface::ReadParameterD(inp, string("T"));
    Restart                  = UserInterface::ReadParameterB(inp, moduleLocation, string("Restrt"));
    tStart                   = UserInterface::ReadParameterI(inp, moduleLocation, string("tStart"));
    tRestartWrite            = UserInterface::ReadParameterI(inp, moduleLocation, string("tRstrt"));
    tFileWrite               = UserInterface::ReadParameterI(inp, moduleLocation, string("FTime"));
    tScreenWrite             = UserInterface::ReadParameterI(inp, moduleLocation, string("STime"));

    inp.close();

    Info::WriteLine();
    Info::WriteLine();

    if (tFileWrite < 1)
    {
        Info::WriteSimple("Illegal value for FWrite; set to 1");
        tFileWrite = 1;
    }
    if (tScreenWrite < 1)
    {
        Info::WriteSimple("Illegal value for STime; set to 1");
        tScreenWrite = 1;
    }

    Eta = iWidth*dx;

    //New
    SD.Initialize();
    SD.ReadInput(InputFileName);

	Cp.ReadInput(*this, InputFileName);
    //end New

	Nphases = Cp.Nphases;
	Ncomp = Cp.Ncomp;

    int ignore1;
    int ignore2;
    struct stat st;
    if(stat(VTKDir.c_str(),&st) != 0)
    {
        ignore1 = system(string("mkdir " + VTKDir).c_str());
        Info::WriteStandard("Create directory", VTKDir);
    }

    struct stat st2;
    if(stat(RawDataDir.c_str(),&st2) != 0)
    {
        ignore2 = system(string("mkdir " + RawDataDir).c_str());
        Info::WriteStandard("Create directory", RawDataDir);
    }

    ignore1 = ignore2; // used to avoid "unused return value" compiler warning
    ignore2 = ignore1;

    Info::WriteLine();
}

Settings& Settings::operator= (const Settings& rhs)
{
    if (this != &rhs) // protect against invalid self-assignment
    {
        thisclassname = rhs.thisclassname;

        //DefaultInputFileName = rhs.DefaultInputFileName;

        SimulationTitle = rhs.SimulationTitle;
        UnitsOfLength = rhs.UnitsOfLength;
        UnitsOfMass = rhs.UnitsOfMass;
        UnitsOfTime = rhs.UnitsOfTime;
        UnitsOfEnergy = rhs.UnitsOfEnergy;

        Nx = rhs.Nx;
        Ny = rhs.Ny;
        Nz = rhs.Nz;

        Nphases = rhs.Nphases;
        Ncomp = rhs.Ncomp;
        dt = rhs.dt;
        nSteps = rhs.nSteps;
        tFileWrite = rhs.tFileWrite;
        tScreenWrite = rhs.tScreenWrite;
        tStart = rhs.tStart;
        Restart = rhs.Restart;
        tRestartWrite = rhs.tRestartWrite;
        nOMP = rhs.nOMP;
        
        omp_set_num_threads(nOMP);
    
        dx = rhs.dx;
        iWidth = rhs.iWidth;
        Eta = rhs.Eta;
    }
    return *this;
}

} //namespace opensim
