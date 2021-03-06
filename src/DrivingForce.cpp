#include "Info.h"
#include "DrivingForce.h"
#include "Tools/UserInterface.h"
#include "VTK.h"
#include "GrainInfo.h"
#include "Settings.h"
#include "PhaseField.h"
#include "InterfaceMobility.h"
#include "InterfaceEnergy.h"
#include "BoundaryConditions.h"
#include <random>
namespace opensim
{
    using namespace std;
    DrivingForce::DrivingForce(const Settings& locSettings, const int boundary)
    {
        this->Initialize(locSettings, boundary);
        this->ReadInput();
    }

    void DrivingForce::Initialize(const Settings& locSettings,
            const int boundary)
    {
        thisclassname = "DrivingForce";
        //DefaultInputFileName = ProjectInputDir + "DrivingForceInput.opi";

        Nx = locSettings.Nx;
        Ny = locSettings.Ny;
        Nz = locSettings.Nz;

        // Setting default values to be used if ReadInput() is not called
        CutOff = 0.95;
        Averaging = false;
        Range = ((locSettings.iWidth + 1)/2)+1;
        PhiThreshold = 1.0/3.0;
        switch(int(locSettings.iWidth))
        {
            case 3:
            {
                PhiThreshold = 1.0/5.0;
                break;
            }
            case 4:
            {
                PhiThreshold = 1.0/4.0;
                break;
            }
            case 5:
            {
                PhiThreshold = 1.0/3.0;
                break;
            }
            default: PhiThreshold = 1.0/3.0;
        }
        // End setting default values

        Raw.Allocate(Nx, Ny, Nz, Range);

        dGabLimitReachedCounter = 0;
        EnDensUnits = locSettings.UnitsOfEnergy + "/" + locSettings.UnitsOfLength + "^3";

        initialized = true;
        Info::WriteStandard(thisclassname, "Initialized");
    }

    void DrivingForce::ReadInput(string InputFileName)
    {
        fstream inp(InputFileName.c_str(), ios::in);

        if (!inp)
        {
            Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
            exit(1);
        };

        Info::WriteBlankLine();

        Info::WriteLineInsert("Driving Force settings");
        Info::WriteStandard("Source", InputFileName.c_str());

        int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

        CutOff = UserInterface::ReadParameterD(inp, moduleLocation, string("CutOff"), false, CutOff);
        Averaging = UserInterface::ReadParameterB(inp, moduleLocation, string("Averg"), false, "Yes");
        Range = UserInterface::ReadParameterI(inp, moduleLocation, string("Range"), false, Range);
        PhiThreshold = UserInterface::ReadParameterD(inp, moduleLocation, string("Thresh"), false, PhiThreshold);

        Info::WriteLine();

        inp.close();
    }

    void DrivingForce::Average(PhaseField& Phase, BoundaryConditions& BC)
    {
        if(Averaging)
        {
            SetBoundaryConditions(BC);
            CollectAverage(Phase);
            SetBoundaryConditions(BC);
            DistributeAverage(Phase);
        }
    }

    void DrivingForce::WriteVTKforPhases(PhaseField& Phi, const int tStep) const
    {
        /**This function will write a readable output in VTS format, summarizing
         * the driving forces of all individual grains as the driving force between
         * individual phase fields.*/

        stringstream converter;
        string phases = converter.str();

        stringstream buffer;
        std::vector<int> DataTypes {PDScalars};

        VTK::WriteHeader(buffer, Nx, Ny, Nz);
        VTK::WriteBeginPointData(buffer, DataTypes);

        for(int alpha = 0; alpha < Phi.Nphases; alpha++)
        for(int beta = alpha+1; beta < Phi.Nphases; beta++)
        {
            buffer << "<DataArray type = \"Float64\" Name = \"dG("
                << alpha << "," << beta << ")\" "
                << "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
            for(int k = 0; k < Nz; ++k)
            for(int j = 0; j < Ny; ++j)
            for(int i = 0; i < Nx; ++i)
            {
                double tempdG = 0.0;
                for(auto it1 = Phi.Fields(i,j,k).cbegin();
                        it1 < Phi.Fields(i,j,k).cend(); ++it1)
                for(auto it2 = Phi.Fields(i,j,k).cbegin();
                        it2 < Phi.Fields(i,j,k).cend(); ++it2)
                if((Phi.FieldsStatistics[it1->index].Phase == alpha)
                and(Phi.FieldsStatistics[it2->index].Phase == beta))
                {
                    tempdG += Raw(i,j,k).get_asym(it1->index, it2->index);
                }
                buffer << tempdG << "\n";
            }
            buffer << "</DataArray>" << endl;
        }
        VTK::WriteEndPointData(buffer);
        VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
        VTK::WriteToFile(buffer, "DrivingForce", tStep);
    }

    void DrivingForce::PrintDiagnostics()
    {
        if (dGabLimitReachedCounter)
        {
            std::string message = "\nThe driving force value has been limited " +
                                std::to_string(static_cast <long int> (dGabLimitReachedCounter)) + " times\n"
                                + "Max value: " + std::to_string(static_cast <long double> (MAXDrivingForce))
                                + " " + EnDensUnits + "\n"
                                + "Max overshoot: " + std::to_string(static_cast <long double> (MAXdGabOvershoot)) + " times\n";
            Info::WriteStandard(thisclassname + "::PrintDiagnostics()", message, false);
            dGabLimitReachedCounter = 0;
            MAXdGabOvershoot        = 0;
            MAXDrivingForce         = 0;
        }
        else
        {
            std::string message = "No Driving Force limiting was needed!";
            Info::WriteSimple(message, false);
        }
    }
}