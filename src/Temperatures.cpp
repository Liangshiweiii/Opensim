#include "Temperatures.h"
#include "Info.h"
#include "Tools/UserInterface.h"
#include "Settings.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "VTK.h"
#include "Velocities.h"

namespace opensim
{
    using namespace std;
    Temperature::Temperature(const Settings& locSettings, const std::string InputFileName)
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

    void Temperature::Initialize(const Settings& locSettings,unsigned int boundary)
    {
        thisclassname = "Temperature";
        //DefaultInputFileName = ProjectInputDir + "TemperatureInput.opi";

        Nx = locSettings.Nx;
        Ny = locSettings.Ny;
        Nz = locSettings.Nz;
        Nphases = locSettings.Nphases;
        dx = locSettings.dx;

        Tx.Allocate(Nx, Ny, Nz, boundary);
        TxDx.Allocate(Nx, Ny, Nz, boundary);
        qdot.Allocate(Nx, Ny, Nz, boundary);

        dT_dr.set_to_zero();
        r0.set_to_zero();
        initialized = true;
        Info::WriteStandard(thisclassname, "Initialized");
    }

    void Temperature::ReadInput(const string InputFileName)
    {
        fstream inp(InputFileName.c_str(), ios::in);

        if (!inp)
        {
            Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
            exit(1);
        };

        Info::WriteBlankLine();
        Info::WriteLineInsert("Temperature properties");
        Info::WriteStandard("Source", InputFileName.c_str());

        int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
        
        T0 = UserInterface::ReadParameterD(inp, moduleLocation, string("T0"));

        r0[0] = UserInterface::ReadParameterD(inp, moduleLocation, string("R0X"), false, 0.0);
        r0[1] = UserInterface::ReadParameterD(inp, moduleLocation, string("R0Y"), false, 0.0);
        r0[2] = UserInterface::ReadParameterD(inp, moduleLocation,  string("R0Z"), false, 0.0);

        dT_dr[0] = UserInterface::ReadParameterD(inp, moduleLocation, string("DT_DRX"), false, 0.0);
        dT_dr[1] = UserInterface::ReadParameterD(inp, moduleLocation, string("DT_DRY"), false, 0.0);
        dT_dr[2] = UserInterface::ReadParameterD(inp, moduleLocation, string("DT_DRZ"), false, 0.0);

        dT_dt = UserInterface::ReadParameterD(inp, moduleLocation, string("DT_Dt"));

        Info::WriteLine();

        inp.close();
    }

    void Temperature::SetBoundaryConditions(const BoundaryConditions& BC)
    {
        BC.SetX(Tx);
        BC.SetY(Tx);
        BC.SetZ(Tx);
    }

    void Temperature::Read(const int tStep)
    {
        string FileName = UserInterface::MakeFileName(RawDataDir,"Temperature_", tStep, ".dat");

        fstream inp(FileName.c_str(), ios::in | ios::binary);

        if (!inp)
        {
            Info::WriteExit("File \"" + FileName + "\" could not be opened", thisclassname, "Read()");
            exit(1);
        };

        inp.read(reinterpret_cast<char*>(&Tx[0]), Tx.tot_size()*sizeof(double));
    }

    void Temperature::SetInitial(const BoundaryConditions& BC, const PhaseField& Phase,
                                                                    const int PhaseIndex)
    {
        RefVolume = Phase.FieldsStatistics[PhaseIndex].Volume;
        SetInitial(BC);
    }

    void Temperature::SetInitial(const BoundaryConditions& BC)
    {
        const long int boundary = Tx.Bcells();
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,boundary,)
        {
            Tx(i,j,k) = T0 + (dT_dr[0]*(i - r0[0]) +
                            dT_dr[1]*(j - r0[1]) +
                            dT_dr[2]*(k - r0[2]))*dx;
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        SetBoundaryConditions(BC);
    }

    void Temperature::Set(const BoundaryConditions& BC, const double dt)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
        {
            Tx(i,j,k) += dt * dT_dt;
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        SetBoundaryConditions(BC);
    }

    void Temperature::Set(const BoundaryConditions& BC, const PhaseField& Phase,
            const double LH, const double Cp, const int PhaseIndex, const double dt)
    {
        double dTtransf = LH*(RefVolume - Phase.FieldsStatistics[PhaseIndex].Volume)/(Cp*(Nx*Ny*Nz));

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
        {
            Tx(i,j,k) += dTtransf + dt * dT_dt;
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        SetBoundaryConditions(BC);

        RefVolume = Phase.FieldsStatistics[PhaseIndex].Volume;
    }

    void Temperature::Set(const BoundaryConditions& BC, const PhaseField& Phase,
            const InterfaceField& Psi,
            const double LH, const double Cp, const int PhaseIndex, const double dt)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,0,)
        if(Phase.Fields(i,j,k).flag)
        {
            double dF_dt = 0.0;
            for(auto it = Phase.FieldsDot(i,j,k).cbegin(); it < Phase.FieldsDot(i,j,k).cend(); ++it)
            {
                if(it->indexA == PhaseIndex)
                {
                    dF_dt += it->value;
                }
                else if(it->indexB == PhaseIndex)
                {
                    dF_dt -= it->value;
                }
            }
            Tx(i,j,k) += dt * (dF_dt*LH/Cp + dT_dt);
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        SetBoundaryConditions(BC);
    }
    void Temperature::Set(const BoundaryConditions& BC, const double dTx,
                                                const double Temperature, const double dt)
    {
        const long int boundary = Tx.Bcells();
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Tx,boundary,)
        {
            Tx(i,j,k) += dTx * dt;
            if (dTx == 0.0) Tx(i,j,k) = Temperature;
            if (dTx > 0.0 and Tx(i,j,k) > Temperature)  Tx(i,j,k) = Temperature;
            if (dTx < 0.0 and Tx(i,j,k) < Temperature)  Tx(i,j,k) = Temperature;
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        SetBoundaryConditions(BC);
    }

    void Temperature::Write(const int tStep) //should be const
    {
        string FileName = UserInterface::MakeFileName(RawDataDir,"Temperature_", tStep, ".dat");

        ofstream out(FileName.c_str(), ios::out | ios::binary);

        if (!out)
        {
            Info::WriteExit("File \"" + FileName + "\" could not be created", thisclassname, "Write()");
            exit(1);
        };

        out.write(reinterpret_cast<char*>(&Tx[0]), Tx.tot_size()*sizeof(double));
    }

    
}