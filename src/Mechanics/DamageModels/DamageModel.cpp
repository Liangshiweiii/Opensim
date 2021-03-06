

#include "Settings.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "Mechanics/DamageModels/DamageModel.h"
#include "Mechanics/Storages/DamageProperties.h"
#include "PhaseField.h"
#include "Tools.h"
#include "Info.h"
#include "Tools/UserInterface.h"
#include "VTK.h"

namespace opensim
{

using namespace std;

DamageModel::DamageModel(const Settings& locSettings, const std::string InputFileName)
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

DamageModel::~DamageModel(void)
{
    fftw_destroy_plan(ForwardPlanKAPPA);
    fftw_destroy_plan(BackwardPlanKAPPANL);

    delete[] damage;
    delete[] kappaF;
    delete[] kappaNL;
    delete[] kappaNLF;

    for(int n = 0; n < 3; n++)
    {
        delete[] Q[n];
    }
}

void DamageModel::ReadInput(const std::string IFileName)
{
    if(IFileName.empty() == true)
    {
        Info::WriteExit("No input file for damage parameters given.", thisclassname, "ReadInput()");
        exit(1);
    }

	Info::WriteLineInsert("Damage Model");
	Info::WriteStandard("Source", IFileName);

    fstream inp(IFileName.c_str(), ios::in);

     if (!inp)
     {
         Info::WriteExit("File \"" + IFileName + "\" could not be opened", thisclassname, "ReadInput()");
         exit(1);
     };
     
     int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
     
     alphaNL = UserInterface::ReadParameterD(inp, moduleLocation, string("alphanonlocal"));
     maxDmg  = UserInterface::ReadParameterD(inp, moduleLocation, string("maxDmg"), false, 0.95);

     for(int n = 0; n < Nphases; n++)
     {
         stringstream converter;
         converter << n;
         string counter = converter.str();
         Info::WriteBlankLine();
         damageflag[n]  = UserInterface::ReadParameterB(inp, moduleLocation, string("damageflag_") + counter);

         k0_linear[n]   = UserInterface::ReadParameterD(inp, moduleLocation, string("k0_linear_") + counter);
         kc_linear[n]   = UserInterface::ReadParameterD(inp, moduleLocation, string("kc_linear_") + counter);

         Info::WriteLine("-");
     }
}

void DamageModel::Initialize(const Settings& locSettings)
{
    //DefaultInputFileName = ProjectInputDir + "DamageInput.opi";
    thisclassname = "DamageModel";

    Nx     = locSettings.Nx;
    Ny     = locSettings.Ny;
    Nz     = locSettings.Nz;
    Nzc    = (locSettings.Nz)/2+1;
    rlSIZE = Nx*Ny*Nz;
    rcSIZE = Nx*Ny*Nzc;
    dx = locSettings.dx;
    iLx = 1.0/(double(Nx)*locSettings.dx);
    iLy = 1.0/(double(Ny)*locSettings.dx);
    iLz = 1.0/(double(Nz)*locSettings.dx);
    Nphases = locSettings.Nphases;

    damageflag.Allocate(Nphases);
    k0_linear.Allocate(Nphases);
    kc_linear.Allocate(Nphases);


    Norm = 1.0/double(rlSIZE);
    for(int n = 0; n < 3; n++)
    {
        Q[n] = new double [rcSIZE] ();
    }
    // Recalculate wave vector
    QXYZ();

    // Damage parameter Initialization
    kappa    = new double[rlSIZE] ();
    kappaNL  = new double[rlSIZE] ();
    damage   = new double[rlSIZE] ();

    for(int i = 0; i < Nx; i++)
    for(int j = 0; j < Ny; j++)
    for(int k = 0; k < Nz; k++)
    {
        damage  [k + Nz*(j + Ny*i)] = 0.0;
        kappa   [k + Nz*(j + Ny*i)] = 0.0;
        kappaNL [k + Nz*(j + Ny*i)] = 0.0;
    }
    kappaF  = new complex<double> [rcSIZE] ();
    kappaNLF  = new complex<double> [rcSIZE] ();
    for(int i = 0; i < Nx; i++)
    for(int j = 0; j < Ny; j++)
    for(int k = 0; k < Nzc; k++)
    {
        int XYZ = k + Nzc*(j + Ny*i);

        kappaF[XYZ] = 0.0;
        kappaNLF[XYZ] = 0.0;
    }

    // Inizialize fft plans
    ForwardPlanKAPPA = fftw_plan_dft_r2c_3d
                    (Nx, Ny, Nz, kappa,
                     reinterpret_cast<fftw_complex*> (kappaF),
                     FFTW_PATIENT);

    BackwardPlanKAPPANL = fftw_plan_dft_c2r_3d
                    (Nx, Ny, Nz,
                     reinterpret_cast<fftw_complex*> (kappaNLF),
                     kappaNL,
                     FFTW_PATIENT);

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void DamageModel::Reinitialize(const long int newNx, const long int newNy,
                         const long int newNz)
{
    Nx     = newNx;
    Ny     = newNy;
    Nz     = newNz;
    Nzc    = Nz/2+1;
    rlSIZE = Nx*Ny*Nz;
    rcSIZE = Nx*Ny*Nzc;
    iLx = 1.0/(double(Nx)*dx);
    iLy = 1.0/(double(Ny)*dx);
    iLz = 1.0/(double(Nz)*dx);

    // Inizialize fft plants
    fftw_destroy_plan(ForwardPlanKAPPA);
    fftw_destroy_plan(BackwardPlanKAPPANL);

    Norm = 1.0/double(rlSIZE);
    for(int n = 0; n < 3; n++)
    {
        delete[] Q[n];
    }
    for(int n = 0; n < 3; n++)
    {
        Q[n] = new double [rcSIZE] ();
    }
    // Recalculate wave vector
    QXYZ();

    delete[] kappa;
    delete[] kappaNL;
    delete[] damage;

    kappa    = new double[rlSIZE] ();
    kappaNL  = new double[rlSIZE] ();
    damage   = new double[rlSIZE] ();

    for(int i = 0; i < Nx; i++)
    for(int j = 0; j < Ny; j++)
    for(int k = 0; k < Nz; k++)
    {
        damage  [k + Nz*(j + Ny*i)] = 0.0;
        kappa   [k + Nz*(j + Ny*i)] = 0.0;
        kappaNL [k + Nz*(j + Ny*i)] = 0.0;
    }

    delete[] kappaF;
    delete[] kappaNLF;

    kappaF  = new complex<double> [rcSIZE] ();
    kappaNLF  = new complex<double> [rcSIZE] ();
    for(int i = 0; i < Nx; i++)
    for(int j = 0; j < Ny; j++)
    for(int k = 0; k < Nzc; k++)
    {
        int XYZ = k + Nzc*(j + Ny*i);

        kappaF[XYZ] = 0.0;
        kappaNLF[XYZ] = 0.0;
    }

    ForwardPlanKAPPA = fftw_plan_dft_r2c_3d
                    (Nx, Ny, Nz, kappa,
                     reinterpret_cast<fftw_complex*> (kappaF),
                     FFTW_MEASURE);

    BackwardPlanKAPPANL = fftw_plan_dft_c2r_3d
                    (Nx, Ny, Nz,
                     reinterpret_cast<fftw_complex*> (kappaNLF),
                     kappaNL,
                     FFTW_MEASURE);

    Info::WriteStandard(thisclassname, "Reinitialized");
}

void DamageModel::Solve(ElasticProperties& EP, DamageProperties& DP,PhaseField& Phi, BoundaryConditions& BC)
{
    // Plastic equivalent strain
    F_SetEquivalentPlasticStrain(DP);
    ExecuteForwardTRAFO_KAPPA();
    F_CalculateNLEquationPerlings(EP);
    // Nonlocal damage to real space
    ExecuteBackwardTRAFO_KAPPANL();
    // Normalize KAPPANL after backward transformation
    F_NormalizeKAPPANL();
    // Calculate damage surface
    F_EValuateDamageSurfacePerlingsDuctile2Intf(Phi,EP,DP);           // Change DmgEvo here

    // Couple brittle and ductile damage
    F_CombineBrittleDuctileDamage(DP);

    // Compare damage values
    DP.CompareDamage();
    DP.SaveDamageN();
    DP.SetBoundaryConditions(BC);
}

void DamageModel::F_CalculateNLEquationPerlings(ElasticProperties& EP)
{
    OMP_PARALLEL_FOR_LOOP_BEGIN(XYZ,0,rcSIZE,)
    {
        double QAMOUNT = Q[0][XYZ]*Q[0][XYZ] + Q[1][XYZ]*Q[1][XYZ] + Q[2][XYZ]*Q[2][XYZ];
        kappaNLF[XYZ] = kappaF[XYZ]/(1.0 + alphaNL*QAMOUNT);
    }
    OMP_PARALLEL_FOR_LOOP_END
}


void DamageModel::QXYZ(void)
{
    // Calculation of frequency vector in fftw vector scheme (considering symmetries)
    for(int i = 0; i < Nx ; i++)
    for(int j = 0; j < Ny ; j++)
    for(int k = 0; k < Nzc; k++)
    {
        int XYZ = k + Nzc*(j + Ny*i);

        Q[0][XYZ] = 2.0*Pi*iLx*(double(i)*(i <= Nx/2) - (Nx-double(i))*(i > Nx/2));
        Q[1][XYZ] = 2.0*Pi*iLy*(double(j)*(j <= Ny/2) - (Ny-double(j))*(j > Ny/2));
        Q[2][XYZ] = 2.0*Pi*iLz*(double(k)*(k <= Nz/2) - (Nz-double(k))*(k > Nz/2));
    }
}

void DamageModel::ExecuteForwardTRAFO_KAPPA(void)
{
    fftw_execute(ForwardPlanKAPPA);
}

void DamageModel::ExecuteBackwardTRAFO_KAPPANL(void)
{
    fftw_execute(BackwardPlanKAPPANL);
}

void DamageModel::F_EValuateDamageSurfacePerlingsDuctile2Intf(PhaseField& Phi, ElasticProperties& EP,DamageProperties& DP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DP.EffectiveDamage_d,0,)
    {
        DP.EffectiveDamage_d(i,j,k) = 0;
        if (Phi.Interface(i,j,k))
        {
            for(auto alpha = Phi.Fields(i, j, k).cbegin(); 
                     alpha != Phi.Fields(i, j, k).cend(); ++alpha)
            {
                if (alpha->value != 0)
                {
                    int pIndex = alpha->index;
                    int thPhaseIndex = Phi.FieldsStatistics[pIndex].Phase;
                    if (damageflag[thPhaseIndex])
                    {
                        int XYZ = k + Nz*(j + Ny*i);
                        double k0 = k0_linear[thPhaseIndex];
                        double kc = kc_linear[thPhaseIndex];
                        double kk = kappaNL[XYZ];

                        double DC = maxDmg;

                        if (kk < k0)
                        {
                            damage[XYZ] = 0.0;
                        }
                        else if (kk < kc)
                        {
                            damage[XYZ] = DC*(kk-k0)/(kc-k0);
                        }
                        else
                        {
                            damage[XYZ] = DC;
                        }
                        DP.EffectiveDamage_d(i,j,k) += alpha->value*damage[XYZ];
                    }
                }
            } // end alpha
        }
        else
        {
            int locIndex = Phi.Fields(i,j,k).front().index;
            int thPhaseIndex = Phi.FieldsStatistics[locIndex].Phase;
            if (damageflag[thPhaseIndex])
            {
                int XYZ = k + Nz*(j + Ny*i);
                double k0 = k0_linear[thPhaseIndex];
                double kc = kc_linear[thPhaseIndex];
                double kk = kappaNL[XYZ];

                double DC = maxDmg;

                if (kk < k0)
                {
                    damage[XYZ] = 0.0;
                }
                else if (kk < kc)
                {
                    damage[XYZ] = DC*(kk-k0)/(kc-k0);
                }
                else
                {
                    damage[XYZ] = DC;
                }
                DP.EffectiveDamage_d(i,j,k) = damage[XYZ];
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DamageModel::F_CombineBrittleDuctileDamage(DamageProperties& DP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DP.EffectiveDamage,0,)
    {
        DP.EffectiveDamage(i,j,k) = DP.EffectiveDamage_b(i,j,k) + DP.EffectiveDamage_d(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DamageModel::F_SetEquivalentPlasticStrain(DamageProperties& DP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DP.peeq,0,)
    {
        int XYZ = k + Nz*(j + Ny*i);
        kappa[XYZ] = DP.peeq(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DamageModel::F_NormalizeKAPPANL()
{
    OMP_PARALLEL_FOR_LOOP_BEGIN(XYZ,0,rlSIZE,)
    {
        kappaNL[XYZ] *= Norm;
    }
    OMP_PARALLEL_FOR_LOOP_END
}

void DamageModel::WriteVTK(int tStep)
{
    stringstream outbufer;

    outbufer << "# vtk DataFile Version 3.0\n";
    outbufer << "Kappa\n";
    outbufer << "ASCII\n";
    outbufer << "DATASET STRUCTURED_GRID\n";
    outbufer << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    outbufer << "POINTS " <<  Nx*Ny*Nz << " double\n";

    for(double k = 0; k < Nz; ++k)
    for(double j = 0; j < Ny; ++j)
    for(double i = 0; i < Nx; ++i)
    {
        outbufer << i << " " << j << " " << k << "\n";
    }
    outbufer << "\n";
    outbufer << "POINT_DATA " << Nx*Ny*Nz << "\n";

    outbufer << "SCALARS kappa" << " double 1\n";
    outbufer << "LOOKUP_TABLE default\n";

    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        int XYZ = k + Nz*(j + Ny*i);
        outbufer << kappa[XYZ] <<"\n";
    }

    outbufer << "SCALARS kappaNL" << " double 1\n";
    outbufer << "LOOKUP_TABLE default\n";

    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        int XYZ = k + Nz*(j + Ny*i);
        outbufer << kappaNL[XYZ] <<"\n";
    }

    string FileName = UserInterface::MakeFileName(VTKDir,"DamageModel_", tStep, ".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbufer.rdbuf();
    vtk_file.close();
}

}// namespace openphase
