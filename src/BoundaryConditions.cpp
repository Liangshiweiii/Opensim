#include "BoundaryConditions.h"
#include "Info.h"
#include "Settings.h"
#include "Tools/UserInterface.h"

namespace opensim
{
    using namespace std;
    BoundaryConditions::BoundaryConditions(const Settings& locSettings, const std::string InputFileName)
    {
        this->Initialize(locSettings);

        if(InputFileName == "NONE")
        {
            this->ReadInput();
        }
        else
        {
            this->ReadInput(InputFileName);
        }
    }

    void BoundaryConditions::Initialize(const Settings& locSettings)
    {
        thisclassname = "BoundaryConditions";
        //DefaultInputFileName = ProjectInputDir + "BoundaryConditions.opi";

        BC0X = Periodic;
        BCNX = Periodic;
        BC0Y = Periodic;
        BCNY = Periodic;
        BC0Z = Periodic;
        BCNZ = Periodic;
        initialized = true;
        Info::WriteStandard(thisclassname, "Initialized");
    }

    void BoundaryConditions::ReadInput(const string InputFileName)
    {
        Info::WriteBlankLine();
        Info::WriteLineInsert("Boundary conditions");
        Info::WriteStandard("Source", InputFileName);

        fstream inp(InputFileName.c_str(), ios::in);

        if (!inp)
        {
            Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
            exit(1);
        };
        int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
        
        string BC0Xstring = UserInterface::ReadParameterF(inp, moduleLocation, string("BC0X"));
        string BCNXstring = UserInterface::ReadParameterF(inp, moduleLocation, string("BCNX"));
        string BC0Ystring = UserInterface::ReadParameterF(inp, moduleLocation, string("BC0Y"));
        string BCNYstring = UserInterface::ReadParameterF(inp, moduleLocation, string("BCNY"));
        string BC0Zstring = UserInterface::ReadParameterF(inp, moduleLocation, string("BC0Z"));
        string BCNZstring = UserInterface::ReadParameterF(inp, moduleLocation, string("BCNZ"));

        inp.close();

        std::vector<std::string> BCTypes;
        BCTypes.push_back("Periodic");
        BCTypes.push_back("NoFlux");
        BCTypes.push_back("Free");
        BCTypes.push_back("Fixed");

        int check = 0;

        for (unsigned int i = 0; i < BCTypes.size(); i++)
        {
            if(BC0Xstring == BCTypes[i])
            {
                BC0X = i;
                check ++;
            }
            if(BCNXstring == BCTypes[i])
            {
                BCNX = i;
                check ++;
            }
            if(BC0Ystring == BCTypes[i])
            {
                BC0Y = i;
                check ++;
            }
            if(BCNYstring == BCTypes[i])
            {
                BCNY = i;
                check ++;
            }
            if(BC0Zstring == BCTypes[i])
            {
                BC0Z = i;
                check ++;
            }
            if(BCNZstring == BCTypes[i])
            {
                BCNZ = i;
                check ++;
            }
        }

        if (check != 6)
        {
            Info::WriteExit(string("Illegal Boundary Conditions Input. ") +
                            string("Legal boundary conditions are ") +
                            string("Periodic, NoFlux, Free or Fixed."),
                            thisclassname, "ReadInput()");
            exit(1);
        }

        Info::WriteLine();
    }
}