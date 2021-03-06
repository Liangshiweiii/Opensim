#include "Diffusion.h"
#include "Settings.h"
#include "Info.h"
#include "GrainInfo.h"
#include "DrivingForce.h"
#include "PhaseField.h"
#include "Compositions.h"
#include "Temperatures.h"
#include "InterfaceMobility.h"
#include "PhysicalConstants.h"
#include "BoundaryConditions.h"
#include "Tools/UserInterface.h"

namespace opensim
{
    using namespace std;
    Diffusion::Diffusion(
            const Settings& locSettings, const int boundary)
    {
        this->Initialize(locSettings, boundary);
        this->ReadInput();

        
    }
    void Diffusion::Initialize(const Settings& locSettings,
            const int boundary)
    {
        thisclassname = "Diffsion";
        //DefaultInputFileName = ProjectInputDir + "DiffusionInput.opi";

        Nphases = locSettings.Nphases;
        dx      = locSettings.dx;
        dx_1    = 1.0/dx;
        dx_2    = 1.0/(dx*dx);
        Eta     = locSettings.Eta;
        Nx = locSettings.Nx;
        Ny = locSettings.Ny;
        Nz = locSettings.Nz;
        Comp = 0;
        R = PhysicalConstants::R;
        TotalMass = 0.0;

        Precision = DBL_EPSILON;

        dMu.Allocate(Nx, Ny, Nz, {Nphases}, boundary);
        DC.Allocate(Nx, Ny, Nz, {Nphases}, boundary);
        PhaseFractions.Allocate(Nx, Ny, Nz, {Nphases}, boundary);

        Pt.Allocate({Nphases, Nphases});
        Ts.Allocate({Nphases, Nphases});
        Cs.Allocate({Nphases, Nphases});
        mL.Allocate({Nphases, Nphases});
        mL_1.Allocate({Nphases, Nphases});

        DC0.resize(Nphases);
        AE.resize(Nphases);
        Stoichiometric.resize(Nphases);
        S.resize(Nphases);

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,dMu,0,)
        {
            for (int n = 0; n < Nphases; ++n)
            {
                dMu(i,j,k)({n}) = 0.0;
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        LimitingNeeded = false;

        DStencil.Set(LaplacianStencil27, dx);

        initialized = true;
        Info::WriteStandard(thisclassname, "Initialized");
    }

    void Diffusion::ReadInput(const string InputFileName)
    {
        fstream inp(InputFileName.c_str(), ios::in);

        if (!inp)
        {
            Info::WriteExit("File " + InputFileName + " could not be opened",thisclassname, "ReadInput()");
            exit(1);
        };

        Info::WriteLineInsert("Diffsion properties");
        Info::WriteStandard("Source", InputFileName);
        
        int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
        
        for(int n = 0; n < Nphases; ++n)
        {
            stringstream converter;
            converter << n;
            string counter = converter.str();
            cout << endl;
            Stoichiometric[n] = UserInterface::ReadParameterB(inp, moduleLocation, string("Flag_") + counter);
            DC0[n]      = UserInterface::ReadParameterD(inp, moduleLocation, string("DC_") + counter);
            AE[n]       = UserInterface::ReadParameterD(inp, moduleLocation, string("AE_") + counter);
            S[n]        = UserInterface::ReadParameterD(inp, moduleLocation, string("EF_") + counter);
        }
        ThereAreStoichiometricPhases = false;
        for(int n = 0; n < Nphases; ++n)
        {
            if (Stoichiometric[n])
            {
                ThereAreStoichiometricPhases = true;
                break;
            }
        }

        for(int n =   0; n < Nphases-1; ++n)
        for(int m = n+1; m < Nphases  ; ++m)
        {
            stringstream converter;
            converter << n << "_" << m;
            string counter = converter.str();
            Cs({n, m})    = Cs({m, n})   = UserInterface::ReadParameterD(inp, moduleLocation, string("Cs_") + counter);
            Ts({n, m})    = Ts({m, n})   = UserInterface::ReadParameterD(inp, moduleLocation, string("Ts_") + counter);
            mL({n, m})    = UserInterface::ReadParameterD(inp, moduleLocation, string("ML_") + counter);
            string tmp  = counter;
            counter[0]  = tmp[2];
            counter[2]  = tmp[0];
            mL({m, n})    = UserInterface::ReadParameterD(inp, moduleLocation, string("ML_") + counter);
        }

        for(int n = 0; n < Nphases; ++n)
        {
            Pt({n, n}) = 1.0;
            mL({n, n}) = 1.0;
            mL_1({n, n}) = 1.0;
            Ts({n, n}) = 0.0;
            Cs({n, n}) = 0.0;
        };

        for(int n = 0; n < Nphases; ++n)
        for(int m = 0; m < Nphases; ++m)
        if(n != m)
        {
            if(mL({n, m}) == 0.0)
            {
                Pt({m, n}) = 0.0;
                mL_1({n, m}) = 0.0;
            }
            else
            {
                mL_1({n, m}) = 1.0/mL({n, m});
            }
            Pt({n, m}) = fabs(mL({m, n}) * mL_1({n, m}));
        }

        inp.close();
        Info::WriteLine();
    }

    void Diffusion::SetPhaseFractions(PhaseField& Phase)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PhaseFractions,PhaseFractions.Bcells(),)
        {
            PhaseFractions(i,j,k) = Phase.Fractions(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    inline void Diffusion::CalculateLocalPhaseConcentrations(PhaseField& Phase,
                                                                Temperature& Tx, Composition& Elements,
                                                                int i, int j, int k)
    {
        int locNphases = 0;
        int locPindex = 0;
        for(int n = 0; n < Nphases; n++)
        if(PhaseFractions(i,j,k)({n}) > 0.0)
        {
            locNphases += 1;
            locPindex = n;
        }

        if(locNphases == 1)
        //if(!Phase.Interface(i,j,k))
        {
            int pIndex = locPindex;
            Elements.Phase(i,j,k)({pIndex, Comp}) = Elements.Total(i,j,k)({Comp});
            if(!Stoichiometric[pIndex])
            {
                for(int n = 0; n < Nphases; ++n)
                if(n != pIndex)
                {
                    Elements.Phase(i,j,k)({n, Comp}) = Cs({n, pIndex}) +
                            (Elements.Phase(i,j,k)({pIndex, Comp}) - Cs({n, pIndex})) *
                            mL({pIndex, n})*mL_1({n, pIndex});
                }
            }
            else
            {
                for(int n = 0; n < Nphases; ++n)
                if(n != pIndex)
                {
                    Elements.Phase(i,j,k)({n, Comp}) = Cs({n, pIndex}) +
                                        (Tx(i,j,k) - Ts({n, pIndex}))*mL_1({n, pIndex});
                }
            }
        }
        else
        {
            Tensor <double, 2> eqCx({Nphases, Nphases});

            for(int n = 0; n < Nphases; ++n)
            for(int m = 0; m < Nphases; ++m)
            if(n != m)
            {
                eqCx({n, m}) = Cs({n, m}) + (Tx(i,j,k) - Ts({n, m}))*mL_1({n, m});
            }

            double dCx = Elements.Total(i,j,k)({Comp});
            for(int n = 0; n < Nphases; ++n)
            {
                Elements.Phase(i,j,k)({n, Comp}) = 0.0;
                eqCx({n,n}) = 0.0;
                double div = 0.0;
                for(int m = 0; m < Nphases; ++m)
                if(n != m)
                {
                    eqCx({n,n}) += eqCx({n, m})*PhaseFractions(i,j,k)({m});
                    div += PhaseFractions(i,j,k)({m});
                }

                if(div > Precision)
                {
                    eqCx({n,n}) /= div;
                }
                else
                {
                    eqCx({n,n}) = Elements.Total(i,j,k)({Comp});
                }
                dCx -= eqCx({n,n})*PhaseFractions(i,j,k)({n});
            }
            //if(i == 6 && j == 0 && k == 2)
            //cout << "1: " << dCx << endl;

            for(int n = 0; n < Nphases; ++n)
            if (!Stoichiometric[n])
            {
                double denom = 0.0;
                double nom = 1.0;
                for(int m = 0; m < Nphases; ++m)
                if (!Stoichiometric[m])
                {
                    //if(i == 6 && j == 0 && k == 2)
                    //cout << "Pt({"<<m<<" " <<n << "}) = " << Pt({m, n}) << endl;
                    denom += Pt({m, n})*PhaseFractions(i,j,k)({m});
                    //if(i == 6 && j == 0 && k == 2)
                    //cout << "3: " <<  Pt({m, n}) << "*" << PhaseFractions(i,j,k)({m}) << endl;
                }
                else
                {
                    nom -= PhaseFractions(i,j,k)({m});
                }

                Elements.Phase(i,j,k)({n, Comp}) = eqCx({n,n});

                if(denom > 0.0)
                {
                    Elements.Phase(i,j,k)({n, Comp}) += nom*dCx/denom;
                    //if(i == 6 && j == 0 && k == 2)
                    //    cout << "2: " << nom*dCx/denom  << " " << nom << " " << denom << endl;
                }
            }
            else
            {
                Elements.Phase(i,j,k)({n, Comp}) = eqCx({n,n});
            }
        }
    }

    void Diffsion::CalculatePhaseConcentrations(PhaseField& Phase,
                                        Temperature& Tx, Composition& Elements)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,Elements.Total.Bcells(),)
        {
            CalculateLocalPhaseConcentrations(Phase, Tx, Elements, i, j, k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    void Diffsion::SetDiffusionCoefficients(PhaseField& Phase,
                                                                    Temperature& Tx)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DC,DC.Bcells(),)
        {
            double invRT  = 1.0/(R*Tx(i,j,k));
            for (int n = 0; n < Nphases; ++n)
            {
                DC(i, j, k)({n}) = DC0[n] * exp(-AE[n] * invRT);
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    void Diffsion::GetDrivingForce(PhaseField& Phase,
                                                        Composition& Elements,
                                                        Temperature& Tx,
                                                        DrivingForce& dGab)
    {
        CalculatePhaseConcentrations(Phase, Tx, Elements);
        GetDrivingForceLoop(Phase, Elements, Tx, dGab);
    }

    void Diffsion::GetDrivingForceLoop(PhaseField& Phase,
                                                        Composition& Elements,
                                                        Temperature& Tx,
                                                        DrivingForce& dGab)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,0,)
        {
            if(Phase.Interface(i,j,k))
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                    alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
            for(auto  beta = alpha + 1;
                    beta < Phase.Fields(i,j,k).cend(); ++beta)
            {
                int index1 = Phase.FieldsStatistics[alpha->index].Phase;
                int index2 = Phase.FieldsStatistics[ beta->index].Phase;

                if(index1 != index2)
                {
                    double dG_AB = 0.5*((Ts({index1, index2}) + mL({index1, index2})*
                                    (Elements.Phase(i,j,k)({index1, Comp}) - Cs({index1, index2})) - Tx(i,j,k))*
                                    (1 - Stoichiometric[index1])*(1 + Stoichiometric[index2])
                                        +
                                    (Ts({index2, index1}) + mL({index2, index1})*
                                    (Elements.Phase(i,j,k)({index2, Comp}) - Cs({index2, index1})) - Tx(i,j,k))*
                                    (1 - Stoichiometric[index2])*(1 + Stoichiometric[index1]))*
                                    (S[index2] - S[index1]);

                    dGab.Raw(i,j,k).add_asym(alpha->index,  beta->index,  dG_AB);
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    double Diffsion::GetDrivingForceAlpha(PhaseField& Phase,
                                                            Composition& Elements,
                                                            Temperature& Tx, int alpha,
                                                            int i, int j, int k)
    {
        double dG = 0.0;
        CalculateLocalPhaseConcentrations(Phase, Tx, Elements, i, j, k);

        for(auto  beta = Phase.Fields(i,j,k).cbegin();
                beta < Phase.Fields(i,j,k).cend(); ++beta)
        {
            int index1 = alpha;
            int index2 = Phase.FieldsStatistics[beta->index].Phase;
            if(index1 != index2 and beta->value != 0.0)
            {
                double dG_AB = 0.5*((Ts({index1, index2}) + mL({index1, index2})*
                                (Elements.Phase(i,j,k)({index1, Comp}) - Cs({index1, index2})) - Tx(i,j,k))*
                                (1 - Stoichiometric[index1])*(1 + Stoichiometric[index2])
                                    +
                                (Ts({index2, index1}) + mL({index2, index1})*
                                (Elements.Phase(i,j,k)({index2, Comp}) - Cs({index2, index1})) - Tx(i,j,k))*
                                (1 - Stoichiometric[index2])*(1 + Stoichiometric[index1]))*
                                (S[index2] - S[index1]);

                dG += dG_AB;
            }
        }
        return dG;
    }


    void Diffsion::RestoreStoichiometric(PhaseField& Phase,
                                                            Temperature& Tx,
                                                            Composition& Elements)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,0,)
        {
            if(Phase.Fields(i,j,k).flag == 1)
            {
                int alpha = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
                if (Stoichiometric[alpha] &&
                    Elements.Total(i,j,k)({Comp}) != Elements.Initial({alpha, Comp}))
                {
                    double mass = 0.0;
                    double dC = Elements.Total(i,j,k)({Comp}) - Elements.Initial({alpha, Comp});
                    for(int x = -1; x <= 1; x++)
                    for(int y = -1; y <= 1; y++)
                    for(int z = -1; z <= 1; z++)
                    if(Phase.Interface(i + x, j + y, k + z) &&
                    (int)i+x >= 0 && (int)i+x < Nx &&
                    (int)j+y >= 0 && (int)j+y < Ny &&
                    (int)k+z >= 0 && (int)k+z < Nz)
                    {
                        for(int beta = 0; beta < Nphases; ++beta)
                        if(!Stoichiometric[beta])
                        {
                            mass += PhaseFractions(i + x, j + y, k + z)({beta});
                        }
                    }

                    if(mass > 0.0)
                    {
                        mass = 1.0/(mass);
                        Elements.Total(i, j, k)({Comp}) -= dC;

                        for(int x = -1; x <= 1; x++)
                        for(int y = -1; y <= 1; y++)
                        for(int z = -1; z <= 1; z++)
                        if(Phase.Interface(i + x, j + y, k + z) &&
                        (int)i+x >= 0 && (int)i+x < Nx &&
                        (int)j+y >= 0 && (int)j+y < Ny &&
                        (int)k+z >= 0 && (int)k+z < Nz)
                        {
                            for(int beta = 0; beta < Nphases; ++beta)
                            if(!Stoichiometric[beta])
                            {
                                Elements.Total(i + x, j + y, k + z)({Comp}) += mass *
                                    dC * PhaseFractions(i + x, j + y, k + z)({beta});
                            }
                            CalculateLocalPhaseConcentrations(Phase, Tx, Elements,i,j,k);
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    void Diffsion::RestoreStoichiometricThreadSafe(PhaseField& Phase,
                                                            Temperature& Tx,
                                                            Composition& Elements)
    {
        int Bcells = Elements.Total.Bcells();
        #ifdef _OPENMP
        int locksize = (Nx+2*Bcells)*(Ny+2*Bcells)*(Nz+2*Bcells);
        std::vector<omp_lock_t> writelock(locksize);
        for (int i = 0; i < locksize; ++i)
        {
            omp_init_lock(&writelock[i]);
        }
        #endif
        Storage3D< double, 0 > dC;
        dC.Allocate(Nx, Ny, Nz, Elements.Total.Bcells());

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,0,)
        {
            if(Phase.Fields(i,j,k).flag == 1)
            {
                int alpha = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
                if (Stoichiometric[alpha] &&
                    Elements.Total(i,j,k)({Comp}) != Elements.Initial({alpha, Comp}))
                {
                    dC(i,j,k) = Elements.Total(i,j,k)({Comp}) - Elements.Initial({alpha, Comp});
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,0,)
        {
            if(Phase.Fields(i,j,k).flag == 1)
            {
                int alpha = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
                if (Stoichiometric[alpha] &&
                    Elements.Total(i,j,k)({Comp}) != Elements.Initial({alpha, Comp}))
                {
                    double mass = 0.0;
                    for(int x = -1; x <= 1; x++)
                    for(int y = -1; y <= 1; y++)
                    for(int z = -1; z <= 1; z++)
                    if(Phase.Interface(i + x, j + y, k + z) &&
                    (int)i+x >= -Bcells && (int)i+x < Nx + Bcells &&
                    (int)j+y >= -Bcells && (int)j+y < Ny + Bcells &&
                    (int)k+z >= -Bcells && (int)k+z < Nz + Bcells)
                    {
                        for(int beta = 0; beta < Nphases; ++beta)
                        if(!Stoichiometric[beta])
                        {
                            mass += PhaseFractions(i + x, j + y, k + z)({beta});
                        }
                    }
                    if(mass > 0.0)
                    {
                        mass = 1.0/(mass);
                        #ifdef _OPENMP
                        omp_set_lock(&writelock[((Ny + 2*Bcells)*(i + Bcells) + j + Bcells)*(Nz + 2*Bcells) + k + Bcells]);
                        #endif
                        Elements.Total(i, j, k)({Comp}) -= dC(i,j,k);
                        CalculateLocalPhaseConcentrations(Phase, Tx, Elements,i, j, k);
                        #ifdef _OPENMP
                        omp_unset_lock(&writelock[((Ny + 2*Bcells)*(i + Bcells) + j + Bcells)*(Nz + 2*Bcells) + k + Bcells]);
                        #endif
                        for(int x = -1; x <= 1; x++)
                        for(int y = -1; y <= 1; y++)
                        for(int z = -1; z <= 1; z++)
                        if(Phase.Interface(i + x, j + y, k + z) &&
                        (int)i+x >= Bcells && (int)i+x < Nx + Bcells &&
                        (int)j+y >= Bcells && (int)j+y < Ny + Bcells &&
                        (int)k+z >= Bcells && (int)k+z < Nz + Bcells)
                        {
                            #ifdef _OPENMP
                            omp_set_lock(&writelock[((Ny + 2*Bcells)*(i + Bcells) + j + Bcells)*(Nz + 2*Bcells) + k + Bcells]);
                            #endif
                            for(int beta = 0; beta < Nphases; ++beta)
                            if(!Stoichiometric[beta])
                            {
                                Elements.Total(i + x, j + y, k + z)({Comp}) += mass *
                                    dC(i,j,k) * PhaseFractions(i + x, j + y, k + z)({beta});
                            }
                            CalculateLocalPhaseConcentrations(Phase, Tx, Elements,i + x, j + y, k + z);
                            #ifdef _OPENMP
                            omp_unset_lock(&writelock[((Ny + 2*Bcells)*(i + Bcells) + j + Bcells)*(Nz + 2*Bcells) + k + Bcells]);
                            #endif
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }


    void Diffsion::Remesh(int newNx, int newNy, int newNz)
    {
        DC.Reallocate(newNx, newNy, newNz);
        dMu.Reallocate(newNx, newNy, newNz);
        PhaseFractions.Reallocate(newNx, newNy, newNz);
        Nx = newNx;
        Ny = newNy;
        Nz = newNz;

        Info::WriteStandard(thisclassname, "Remeshed!");
    }

    void Diffsion::CalculateInterfaceMobility(PhaseField& Phase,
                                Composition& Elements, BoundaryConditions& BC, InterfaceMobility& IMobility)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,Elements.Total.Bcells(),)
        {
            if(Phase.Fields(i,j,k).flag)
            {
                IMobility.clear(i,j,k);

                for(auto alpha  = Phase.Fields(i,j,k).cbegin();
                        alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto beta   = alpha + 1;
                        beta  != Phase.Fields(i,j,k).cend(); ++beta)
                {
                    int index1 = Phase.FieldsStatistics[alpha->index].Phase;
                    int index2 = Phase.FieldsStatistics[ beta->index].Phase;

                    double mob = IMobility.MaxIntMobility(index1, index2);
                    if(!Stoichiometric[index1] || !Stoichiometric[index2])
                    {
                        mob = IMobility.IntMobility(index1, index2);

                        double dS = (S[index1] - S[index2]);
                        double dC = (Elements.Phase(i,j,k)({index1,Comp}) -
                                    Elements.Phase(i,j,k)({index2,Comp}));
                        if(alpha->value*beta->value != 0.0 && dS != 0.0 && dC != 0.0)
                        {
                            double mLtmp = 0.0;
                            if(!Stoichiometric[index1])
                            {
                                mLtmp = min(fabs(mL({index1, index2})), fabs(mL({index2, index1})));
                            }
                            else
                            {
                                mLtmp = mL({index2, index1});
                            }

                            double denominator = mLtmp*Eta*dS*dC;

                            if(denominator != 0.0)
                            {
                                mob = fabs(8.0*(DC(i,j,k)({index1}) + DC(i,j,k)({index2}))/denominator);
                            }
                        }
                    }
                    //IMobility.set_raw(i, j, k, alpha->index,  beta->index, mob);
                    IMobility.set(i, j, k, alpha->index,  beta->index, mob);
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        //IMobility.Average(Phase, BC);
    }

    double Diffsion::SetInitialComposition(PhaseField& Phase,
                                                            Composition& Elements)
    {
        double total = 0.0;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total, Elements.Total.Bcells(),reduction(+:total))
        {
            Elements.Phase(i,j,k) = Elements.Initial;
            Elements.Total(i, j, k)({Comp}) = 0.0;
            for(int alpha = 0; alpha != Nphases; ++alpha)
            {
                Elements.Total(i, j, k)({Comp}) +=
                        Elements.Initial({alpha, Comp})*Phase.Fractions(i,j,k)({alpha});

            }
            if((int)i >= 0 && (int)i < Nx && (int)j >= 0 && (int)j < Ny && (int)k >= 0 && (int)k < Nz)
            {
                total += Elements.Total(i, j, k)({Comp});
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        return total;
    }

    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< NEW >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    void Diffsion::CalculateAntitrappingIncrements(PhaseField& Phase,
                                                            Composition& Elements)
    {
        const double koef = Eta/Pi;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,0,)
        {
            if(Phase.Interface(i,j,k))
            for (int x = -1; x <= 1; ++x)
            for (int y = -1; y <= 1; ++y)
            for (int z = -1; z <= 1; ++z)
            if(Phase.Interface(i+x, j+y, k+z) and
            (x != 0 or y != 0 or z != 0))
            {
                Node locPhaseFields = (Phase.Fields(i,j,k).add_exist(Phase.Fields(i+x, j+y, k+z)))*0.5;
                NodeV locNormals = (Phase.Normals(i,j,k) + Phase.Normals(i+x,j+y,k+z));
                if(locPhaseFields.size() > 1)
                for(auto alpha = locPhaseFields.cbegin();
                        alpha != locPhaseFields.cend(); ++alpha)
                {
                    int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                    //if(!Stoichiometric[pIndexA])
                    {
                        for(auto beta = locPhaseFields.cbegin();
                                beta < locPhaseFields.cend(); ++beta)
                        if(alpha != beta)
                        {
                            int pIndexB = Phase.FieldsStatistics[beta->index].Phase;
                            if (pIndexA != pIndexB and (DC0[pIndexA] + DC0[pIndexB]) != 0.0)

                            {
                                double StencilNorm = 1.0/sqrt(x*x + y*y + z*z);
                                double scaleA = 1.0;
                                double scaleB = 1.0;

                                if(Phase.FieldsStatistics[alpha->index].Stage == 1)
                                {
                                    scaleA = (Phase.FieldsStatistics[alpha->index].Volume/Phase.RefVolume);
                                }
                                if(Phase.FieldsStatistics[ beta->index].Stage == 1)
                                {
                                    scaleB = (Phase.FieldsStatistics[beta->index].Volume/Phase.RefVolume);
                                }

                                double NormalXproj = x * locNormals.get_X(alpha->index, beta->index)*0.5;

                                double NormalYproj = y * locNormals.get_Y(alpha->index, beta->index)*0.5;

                                double NormalZproj = z * locNormals.get_Z(alpha->index, beta->index)*0.5;

                                double Projection = (NormalXproj + NormalYproj + NormalZproj)*StencilNorm;

                                double locPsi = 0.5*(Phase.FieldsDot(i  , j  , k  ).get_asym(alpha->index, beta->index) +
                                                    Phase.FieldsDot(i+x, j+y, k+z).get_asym(alpha->index, beta->index));

                                double locDCalpha = 0.5*(DC(i,j,k)({pIndexA}) + DC(i+x, j+y, k+z)({pIndexA}));
                                double locDCbeta  = 0.5*(DC(i,j,k)({pIndexB}) + DC(i+x, j+y, k+z)({pIndexB}));

                                double localDelta = -LaplacianStencil27[x+1][y+1][z+1]*dx_1*
                                                    sqrt(alpha->value*beta->value)*0.5*
                                                    ((Elements.Phase(i  , j  , k  )({pIndexA, Comp}) +
                                                    Elements.Phase(i+x, j+y, k+z)({pIndexA, Comp})) -
                                                    (Elements.Phase(i  , j  , k  )({pIndexB, Comp}) +
                                                    Elements.Phase(i+x, j+y, k+z)({pIndexB, Comp}))) *
                                                    locPsi * Projection * koef *
                                                    scaleA * scaleB *
                                                        (locDCalpha - locDCbeta)/
                                                        (locDCalpha + locDCbeta);

                                if (localDelta < -Precision)
                                {
                                    Elements.PhaseDot(i, j, k)({pIndexA, Comp}).out += localDelta;
                                }
                                if (localDelta > Precision)
                                {
                                    Elements.PhaseDot(i, j, k)({pIndexA, Comp}).in += localDelta;
                                }
                            }
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    void Diffsion::LimitAntitrappingIncrements(PhaseField& Phase,
                                                            Composition& Elements,
                                                            Temperature& Tx)
    {
        const double koef = Eta/Pi;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,0,)
        {
            if(Elements.Limiting(i, j, k))
            {
                if(Phase.Interface(i,j,k))
                {
                    for (int x = -1; x <= 1; ++x)
                    for (int y = -1; y <= 1; ++y)
                    for (int z = -1; z <= 1; ++z)
                    if(Phase.Interface(i+x, j+y, k+z))
                    {
                        Node locPhaseFields = (Phase.Fields(i,j,k).add_exist(Phase.Fields(i+x, j+y, k+z)))*0.5;
                        if(locPhaseFields.size() > 1)
                        for(auto alpha = locPhaseFields.cbegin();
                                alpha != locPhaseFields.cend(); ++alpha)
                        {
                            int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                            //if(!Stoichiometric[pIndexA])
                            for(auto beta = locPhaseFields.cbegin();
                                    beta < locPhaseFields.cend(); ++beta)
                            if(alpha != beta)
                            {
                                int pIndexB = Phase.FieldsStatistics[beta->index].Phase;
                                if (pIndexA != pIndexB and (DC0[pIndexA] + DC0[pIndexB]) != 0.0)

                                {
                                    double StencilNorm = 1.0/sqrt(x*x + y*y + z*z);
                                    double scaleA = 1.0;
                                    double scaleB = 1.0;

                                    if(Phase.FieldsStatistics[alpha->index].Stage == 1)
                                    {
                                        scaleA = (Phase.FieldsStatistics[alpha->index].Volume/Phase.RefVolume);
                                    }
                                    if(Phase.FieldsStatistics[ beta->index].Stage == 1)
                                    {
                                        scaleB = (Phase.FieldsStatistics[beta->index].Volume/Phase.RefVolume);
                                    }

                                    double NormalXproj = x * 0.5*(Phase.Normals(i  , j  , k  ).get_X(alpha->index, beta->index) +
                                                                Phase.Normals(i+x, j+y, k+z).get_X(alpha->index, beta->index));

                                    double NormalYproj = y * 0.5*(Phase.Normals(i  , j  , k  ).get_Y(alpha->index, beta->index) +
                                                                Phase.Normals(i+x, j+y, k+z).get_Y(alpha->index, beta->index));

                                    double NormalZproj = z * 0.5*(Phase.Normals(i  , j  , k  ).get_Z(alpha->index, beta->index) +
                                                                Phase.Normals(i+x, j+y, k+z).get_Z(alpha->index, beta->index));

                                    double Projection = (NormalXproj + NormalYproj + NormalZproj)*StencilNorm;

                                    double locPsi = 0.5*(Phase.FieldsDot(i  , j  , k  ).get_asym(alpha->index, beta->index) +
                                                        Phase.FieldsDot(i+x, j+y, k+z).get_asym(alpha->index, beta->index));

                                    double locDCalpha = 0.5*(DC(i,j,k)({pIndexA}) + DC(i+x, j+y, k+z)({pIndexA}));
                                    double locDCbeta  = 0.5*(DC(i,j,k)({pIndexB}) + DC(i+x, j+y, k+z)({pIndexB}));

                                    double localDelta = -dx*LaplacianStencil27[x+1][y+1][z+1]*
                                                        sqrt(alpha->value*beta->value)*0.5*
                                                        ((Elements.Phase(i  , j  , k  )({pIndexA, Comp}) +
                                                        Elements.Phase(i+x, j+y, k+z)({pIndexA, Comp})) -
                                                        (Elements.Phase(i  , j  , k  )({pIndexB, Comp}) +
                                                        Elements.Phase(i+x, j+y, k+z)({pIndexB, Comp}))) *
                                                        locPsi * Projection * koef *
                                                        scaleA * scaleB *
                                                            (locDCalpha - locDCbeta)/
                                                            (locDCalpha + locDCbeta);

                                    if (Elements.Norm(i, j, k)({pIndexA, Comp}).in <=
                                        Elements.Norm(i+x, j+y, k+z)({pIndexA, Comp}).out and
                                        localDelta > Precision)
                                    {
                                        double localNorm = Elements.Norm(i, j, k)({pIndexA, Comp}).in;
                                        double localIncrement = localDelta*(1.0 - localNorm);
                                        if(localIncrement > Precision)
                                        {
                                            Elements.PhaseDot(  i,   j,   k)({pIndexA, Comp}).in -=
                                                                            localIncrement;
                                            Elements.PhaseDot(i+x, j+y, k+z)({pIndexA, Comp}).out +=
                                                                            localIncrement;
                                        }
                                    }
                                    if (Elements.Norm(i, j, k)({pIndexA, Comp}).out <
                                        Elements.Norm(i+x, j+y, k+z)({pIndexA, Comp}).in and
                                        localDelta < -Precision)
                                    {
                                        double localNorm = Elements.Norm(i, j, k)({pIndexA, Comp}).out;
                                        double localIncrement = localDelta*(1.0 - localNorm);
                                        if(localIncrement < -Precision)
                                        {
                                            Elements.PhaseDot(  i,   j,  k)({pIndexA, Comp}).out -=
                                                                            localIncrement;
                                            Elements.PhaseDot(i+x, j+y, k+z)({pIndexA, Comp}).in +=
                                                                            localIncrement;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    void Diffsion::LimitAntitrappingIncrementsThreadSafe(PhaseField& Phase,
                                                            Composition& Elements,
                                                            Temperature& Tx)
    {
        const double koef = Eta/Pi;
        int Bcells = Elements.Total.Bcells();
        #ifdef _OPENMP
        int locksize = (Nx+2*Bcells)*(Ny+2*Bcells)*(Nz+2*Bcells);
        std::vector<omp_lock_t> writelock(locksize);
        for (int i = 0; i < locksize; ++i)
        {
            omp_init_lock(&writelock[i]);
        }
        #endif

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,0,)
        {
            if(Elements.Limiting(i, j, k))
            {
                if(Phase.Interface(i,j,k))
                {
                    for(int x = -1; x <= 1; ++x)
                    for(int y = -1; y <= 1; ++y)
                    for(int z = -1; z <= 1; ++z)
                    if(x != 0 or y != 0 or z != 0)
                    if(Phase.Interface(i+x, j+y, k+z))
                    {
                        Node locPhaseFields = (Phase.Fields(i,j,k).add_exist(Phase.Fields(i+x, j+y, k+z)))*0.5;
                        NodeV locNormals = (Phase.Normals(i,j,k) + Phase.Normals(i+x, j+y, k+z));
                        if(locPhaseFields.size() > 1)
                        for(auto alpha = locPhaseFields.cbegin();
                                alpha != locPhaseFields.cend(); ++alpha)
                        {
                            int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                            //if(!Stoichiometric[pIndexA])
                            for(auto beta = locPhaseFields.cbegin();
                                    beta < locPhaseFields.cend(); ++beta)
                            if(alpha != beta)
                            {
                                int pIndexB = Phase.FieldsStatistics[beta->index].Phase;
                                if (pIndexA != pIndexB and (DC0[pIndexA] + DC0[pIndexB]) != 0.0)

                                {
                                    double StencilNorm = 1.0/sqrt(x*x + y*y + z*z);
                                    double scaleA = 1.0;
                                    double scaleB = 1.0;

                                    if(Phase.FieldsStatistics[alpha->index].Stage == 1)
                                    {
                                        scaleA = (Phase.FieldsStatistics[alpha->index].Volume/Phase.RefVolume);
                                    }
                                    if(Phase.FieldsStatistics[ beta->index].Stage == 1)
                                    {
                                        scaleB = (Phase.FieldsStatistics[beta->index].Volume/Phase.RefVolume);
                                    }

                                    double NormalXproj = x * locNormals.get_X(alpha->index, beta->index)*0.5;

                                    double NormalYproj = y * locNormals.get_Y(alpha->index, beta->index)*0.5;

                                    double NormalZproj = z * locNormals.get_Z(alpha->index, beta->index)*0.5;

                                    double Projection = (NormalXproj + NormalYproj + NormalZproj)*StencilNorm;

                                    double locPsi = 0.5*(Phase.FieldsDot(i  , j  , k  ).get_asym(alpha->index, beta->index) +
                                                        Phase.FieldsDot(i+x, j+y, k+z).get_asym(alpha->index, beta->index));

                                    double locDCalpha = 0.5*(DC(i,j,k)({pIndexA}) + DC(i+x, j+y, k+z)({pIndexA}));
                                    double locDCbeta  = 0.5*(DC(i,j,k)({pIndexB}) + DC(i+x, j+y, k+z)({pIndexB}));

                                    double localDelta = -LaplacianStencil27[x+1][y+1][z+1]*dx_1*
                                                        sqrt(alpha->value*beta->value)*0.5*
                                                        ((Elements.Phase(i  , j  , k  )({pIndexA, Comp}) +
                                                        Elements.Phase(i+x, j+y, k+z)({pIndexA, Comp})) -
                                                        (Elements.Phase(i  , j  , k  )({pIndexB, Comp}) +
                                                        Elements.Phase(i+x, j+y, k+z)({pIndexB, Comp}))) *
                                                        locPsi * Projection * koef *
                                                        scaleA * scaleB *
                                                            (locDCalpha - locDCbeta)/
                                                            (locDCalpha + locDCbeta);

                                    if (Elements.Norm(i, j, k)({pIndexA, Comp}).in <=
                                        Elements.Norm(i+x, j+y, k+z)({pIndexA, Comp}).out and
                                        localDelta > Precision)
                                    {
                                        double localNorm = Elements.Norm(i, j, k)({pIndexA, Comp}).in;
                                        double localIncrement = localDelta*(1.0 - localNorm);
                                        if(localIncrement > Precision)
                                        {
                                            Elements.PhaseDot(  i,   j,   k)({pIndexA, Comp}).in -= localIncrement;
                                            #ifdef _OPENMP
                                            omp_set_lock(&writelock[((Ny + 2*Bcells)*(i + Bcells) + j + Bcells)*(Nz + 2*Bcells) + k + Bcells]);
                                            #endif
                                            Elements.PhaseDot(i+x, j+y, k+z)({pIndexA, Comp}).out += localIncrement;
                                            #ifdef _OPENMP
                                            omp_unset_lock(&writelock[((Ny + 2*Bcells)*(i + Bcells) + j + Bcells)*(Nz + 2*Bcells) + k + Bcells]);
                                            #endif
                                        }
                                    }
                                    if (Elements.Norm(i, j, k)({pIndexA, Comp}).out <
                                        Elements.Norm(i+x, j+y, k+z)({pIndexA, Comp}).in and
                                        localDelta < -Precision)
                                    {
                                        double localNorm = Elements.Norm(i, j, k)({pIndexA, Comp}).out;
                                        double localIncrement = localDelta*(1.0 - localNorm);
                                        if(localIncrement < -Precision)
                                        {
                                            Elements.PhaseDot(  i,   j,  k)({pIndexA, Comp}).out -= localIncrement;
                                            #ifdef _OPENMP
                                            omp_set_lock(&writelock[((Ny + 2*Bcells)*(i + Bcells) + j + Bcells)*(Nz + 2*Bcells) + k + Bcells]);
                                            #endif
                                            Elements.PhaseDot(i+x, j+y, k+z)({pIndexA, Comp}).in += localIncrement;
                                            #ifdef _OPENMP
                                            omp_unset_lock(&writelock[((Ny + 2*Bcells)*(i + Bcells) + j + Bcells)*(Nz + 2*Bcells) + k + Bcells]);
                                            #endif
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    void Diffsion::CalculateDiffusionIncrements(PhaseField& Phase,
                                                            Composition& Elements,
                                                            Temperature& Tx)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,0,)
        {
            for(int alpha = 0; alpha != Nphases; ++alpha)
            {
                Elements.PhaseDot(i, j, k)({alpha, Comp}).in = 0.0;
                Elements.PhaseDot(i, j, k)({alpha, Comp}).out = 0.0;
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,0,)
        {
            if(Phase.Fields(i,j,k).flag)
            {
                for(int alpha = 0; alpha != Nphases; ++alpha)
                if(!Stoichiometric[alpha] and PhaseFractions(i, j, k)({alpha}) > 0.0)
                {
                    double invRT  = 1.0e4*Elements.MolarVolume({alpha, Comp})/(R*Tx(i,j,k));
                    for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
                    {
                        int di = ds->di;
                        int dj = ds->dj;
                        int dk = ds->dk;
                        if (PhaseFractions(i+di, j+dj, k+dk)({alpha}) > 0.0)
                        {
                            double invRTxyz  = 1.0e4*Elements.MolarVolume({alpha, Comp})/(R*Tx(i+di, j+dj, k+dk));
                            double localDelta = -ds->weight*
                                    sqrt(PhaseFractions(i, j, k)({alpha}) *
                                        PhaseFractions(i+di, j+dj, k+dk)({alpha}))*

                                0.5*((DC(i,j,k)({alpha}) + DC(i+di, j+dj, k+dk)({alpha}))*    // D * Grad ( C {alpha})
                                        (Elements.Phase(  i,   j,   k)({alpha, Comp}) -
                                        Elements.Phase(i+di, j+dj, k+dk)({alpha, Comp})) +

                                        (invRT*DC(i,j,k)({alpha}) +
                                        invRTxyz*DC(i+di, j+dj, k+dk)({alpha}))*
                                        (dMu(  i,   j,   k)({alpha}) -
                                        dMu(i+di, j+dj, k+dk)({alpha})));          // dMu = 0;

                            if (localDelta < -Precision)
                            {
                                Elements.PhaseDot(i, j, k)({alpha, Comp}).out += localDelta;
                            }
                            if (localDelta > Precision)
                            {
                                Elements.PhaseDot(i, j, k)({alpha, Comp}).in += localDelta;
                            }
                        }
                    }
                }
            }
            else
            {
                int alpha = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
                if(!Stoichiometric[alpha])
                {
                    double invRT  = 1.0e4*Elements.MolarVolume({alpha, Comp})/(R*Tx(i,j,k));
                    for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
                    {
                        int di = ds->di;
                        int dj = ds->dj;
                        int dk = ds->dk;
                        double invRTxyz  = 1.0e4*Elements.MolarVolume({alpha, Comp})/(R*Tx(i+di, j+dj, k+dk));
                        double localDelta = - ds->weight*

                            0.5*((DC(i,j,k)({alpha}) + DC(i+di, j+dj, k+dk)({alpha}))*
                                    (Elements.Phase(  i,   j,   k)({alpha, Comp}) -
                                    Elements.Phase(i+di, j+dj, k+dk)({alpha, Comp})) +

                                    (invRT*DC(i,j,k)({alpha}) +
                                    invRTxyz*DC(i+di, j+dj, k+dk)({alpha}))*
                                    (dMu(  i,   j,   k)({alpha}) -
                                    dMu(i+di, j+dj, k+dk)({alpha})));

                        if (localDelta < -Precision)
                        {
                            Elements.PhaseDot(i, j, k)({alpha, Comp}).out += localDelta;
                        }
                        if (localDelta > Precision)
                        {
                            Elements.PhaseDot(i, j, k)({alpha, Comp}).in += localDelta;
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    void Diffsion::CalculateLimits(PhaseField& Phase,
                                                Composition& Elements, double dt)
    {
        bool ompLimitingNeeded = false;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,Elements.Total.Bcells(),)
        {
            for(int alpha = 0; alpha != Nphases; ++alpha)
            {
                Elements.Norm(i,j,k)({alpha, Comp}).in = 1.0;
                Elements.Norm(i,j,k)({alpha, Comp}).out = 1.0;
            }
            Elements.Limiting(i,j,k) = 0;
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,0,reduction(||:ompLimitingNeeded))
        {
            if(Phase.Fields(i,j,k).flag)
            {
                for(int alpha = 0; alpha != Nphases; ++alpha)
                if(!Stoichiometric[alpha] and PhaseFractions(i,j,k)({alpha}) > 0.0)
                {
                    double dPhaseIn  = Elements.PhaseDot(i,j,k)({alpha, Comp}).in*dt;
                    double dPhaseOut = Elements.PhaseDot(i,j,k)({alpha, Comp}).out*dt;

                    if(dPhaseIn <= Precision) dPhaseIn = 0.0;
                    if(dPhaseOut >= -Precision) dPhaseOut = 0.0;

                    double PhaseOld = Elements.Phase(i,j,k)({alpha, Comp});

                    if(PhaseOld + dPhaseIn > Elements.Max({alpha, Comp}) && dPhaseIn != 0.0)
                    {
                        Elements.Norm(i,j,k)({alpha, Comp}).in =
                                            (Elements.Max({alpha, Comp}) - PhaseOld)/dPhaseIn;
                        ompLimitingNeeded = true;

                        for (int x = -1; x <= 1; ++x)
                        for (int y = -1; y <= 1; ++y)
                        for (int z = -1; z <= 1; ++z)
                        {
                            Elements.Limiting(i+x, j+y, k+z) = 1;
                        }
                    }
                    if(PhaseOld + dPhaseOut < Elements.Min({alpha, Comp}) && dPhaseOut != 0.0)
                    {
                        Elements.Norm(i,j,k)({alpha, Comp}).out =
                                            (Elements.Min({alpha, Comp}) - PhaseOld)/dPhaseOut;
                        ompLimitingNeeded = true;
                        for (int x = -1; x <= 1; ++x)
                        for (int y = -1; y <= 1; ++y)
                        for (int z = -1; z <= 1; ++z)
                        {
                            Elements.Limiting(i+x, j+y, k+z) = 1;
                        }
                    }
                }
            }
            else
            {
                int alpha = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
                if(!Stoichiometric[alpha])
                {
                    double dPhaseIn  = Elements.PhaseDot(i,j,k)({alpha, Comp}).in*dt;
                    double dPhaseOut = Elements.PhaseDot(i,j,k)({alpha, Comp}).out*dt;

                    if(dPhaseIn <= Precision) dPhaseIn = 0.0;
                    if(dPhaseOut >= -Precision) dPhaseOut = 0.0;

                    double PhaseOld = Elements.Phase(i,j,k)({alpha, Comp});

                    if(PhaseOld + dPhaseIn > Elements.Max({alpha, Comp}) && dPhaseIn != 0.0)
                    {
                        Elements.Norm(i,j,k)({alpha, Comp}).in =
                                            (Elements.Max({alpha, Comp}) - PhaseOld)/dPhaseIn;
                        ompLimitingNeeded = true;
                        for (int x = -1; x <= 1; ++x)
                        for (int y = -1; y <= 1; ++y)
                        for (int z = -1; z <= 1; ++z)
                        {
                            Elements.Limiting(i+x, j+y, k+z) = 1;
                        }
                    }
                    if(PhaseOld + dPhaseOut < Elements.Min({alpha, Comp}) && dPhaseOut != 0.0)
                    {
                        Elements.Norm(i,j,k)({alpha, Comp}).out =
                                            (Elements.Min({alpha, Comp}) - PhaseOld)/dPhaseOut;
                        ompLimitingNeeded = true;
                        for (int x = -1; x <= 1; ++x)
                        for (int y = -1; y <= 1; ++y)
                        for (int z = -1; z <= 1; ++z)
                        {
                            Elements.Limiting(i+x, j+y, k+z) = 1;
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        LimitingNeeded = ompLimitingNeeded;
    }

    void Diffsion::LimitDiffusionIncrements(PhaseField& Phase,
                                                            Composition& Elements,
                                                            Temperature& Tx)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,0,)
        {
            if(Elements.Limiting(i, j, k))
            {
                if(Phase.Fields(i,j,k).flag)
                {
                    for(int alpha = 0; alpha != Nphases; ++alpha)
                    if(!Stoichiometric[alpha] and PhaseFractions(i, j, k)({alpha}) > 0.0)
                    {
                        double invRT = 1.0e4*Elements.MolarVolume({alpha, Comp})/(R*Tx(i, j, k));
                        for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
                        {
                            int di = ds->di;
                            int dj = ds->dj;
                            int dk = ds->dk;
                            if (PhaseFractions(i+di, j+dj, k+dk)({alpha}) > 0.0)
                            {
                                double invRTxyz  = 1.0e4*Elements.MolarVolume({alpha, Comp})/(R*Tx(i+di, j+dj, k+dk));
                                double localDelta = -ds->weight*
                                        sqrt(PhaseFractions(i, j, k)({alpha}) *
                                            PhaseFractions(i+di, j+dj, k+dk)({alpha}))*

                                    0.5*((DC(i,j,k)({alpha}) + DC(i+di, j+dj, k+dk)({alpha}))*
                                            (Elements.Phase(  i,   j,   k)({alpha, Comp}) -
                                            Elements.Phase(i+di, j+dj, k+dk)({alpha, Comp})) +

                                            (invRT*DC(i,j,k)({alpha}) +
                                            invRTxyz*DC(i+di, j+dj, k+dk)({alpha}))*
                                            (dMu(  i,   j,   k)({alpha}) -
                                            dMu(i+di, j+dj, k+dk)({alpha})));

                                if (Elements.Norm(i, j, k)({alpha, Comp}).in <=
                                    Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).out)
                                {
                                    double localNorm = Elements.Norm(i, j, k)({alpha, Comp}).in;
                                    double localIncrement = localDelta*(1.0 - localNorm);
                                    if(localIncrement > Precision)
                                    {
                                        Elements.PhaseDot(  i,   j,   k)({alpha, Comp}).in -=
                                                                        localIncrement;
                                    }
                                }
                                if (Elements.Norm(i, j, k)({alpha, Comp}).in >
                                    Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).out)
                                {
                                    double localNorm = Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).out;
                                    double localIncrement = localDelta*(1.0 - localNorm);
                                    if(localIncrement > Precision)
                                    {
                                        Elements.PhaseDot(  i,   j,   k)({alpha, Comp}).in -=
                                                                        localIncrement;
                                    }
                                }

                                if (Elements.Norm(i, j, k)({alpha, Comp}).out <
                                    Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).in)
                                {
                                    double localNorm = Elements.Norm(i, j, k)({alpha, Comp}).out;
                                    double localIncrement = localDelta*(1.0 - localNorm);
                                    if(localIncrement < -Precision)
                                    {
                                        Elements.PhaseDot(  i,   j,   k)({alpha, Comp}).out -=
                                                                        localIncrement;
                                    }
                                }
                                if (Elements.Norm(i, j, k)({alpha, Comp}).out >=
                                    Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).in)
                                {
                                    double localNorm = Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).in;
                                    double localIncrement = localDelta*(1.0 - localNorm);
                                    if(localIncrement < -Precision)
                                    {
                                        Elements.PhaseDot(  i,   j,   k)({alpha, Comp}).out -=
                                                                        localIncrement;
                                    }

                                }
                            }
                        }
                    }
                }
                else
                {
                    int alpha = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
                    if(!Stoichiometric[alpha])
                    {
                        double invRT  = 1.0e4*Elements.MolarVolume({alpha, Comp})/(R*Tx(i, j, k));
                        for(auto ds = DStencil.begin(); ds != DStencil.end(); ds++)
                        {
                            int di = ds->di;
                            int dj = ds->dj;
                            int dk = ds->dk;
                            double invRTxyz  = 1.0e4*Elements.MolarVolume({alpha, Comp})/(R*Tx(i+di, j+dj, k+dk));
                            double localDelta = - ds->weight*

                                0.5*((DC(i,j,k)({alpha}) + DC(i+di, j+dj, k+dk)({alpha}))*
                                        (Elements.Phase(  i,   j,   k)({alpha, Comp}) -
                                        Elements.Phase(i+di, j+dj, k+dk)({alpha, Comp})) +

                                        (invRT*DC(i,j,k)({alpha}) +
                                        invRTxyz*DC(i+di, j+dj, k+dk)({alpha}))*
                                        (dMu(  i,   j,   k)({alpha}) -
                                        dMu(i+di, j+dj, k+dk)({alpha})));

                            if (Elements.Norm(i, j, k)({alpha, Comp}).in <=
                                Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).out)
                            {
                                double localNorm = Elements.Norm(i, j, k)({alpha, Comp}).in;
                                double localIncrement = localDelta*(1.0 - localNorm);
                                if(localIncrement > Precision)
                                {
                                    Elements.PhaseDot(  i,   j,   k)({alpha, Comp}).in -=
                                                                    localIncrement;
                                }
                            }
                            if (Elements.Norm(i, j, k)({alpha, Comp}).in >
                                Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).out)
                            {
                                double localNorm = Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).out;
                                double localIncrement = localDelta*(1.0 - localNorm);
                                if(localIncrement > Precision)
                                {
                                    Elements.PhaseDot(  i,   j,   k)({alpha, Comp}).in -=
                                                                    localIncrement;
                                }
                            }

                            if (Elements.Norm(i, j, k)({alpha, Comp}).out <
                                Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).in)
                            {
                                double localNorm = Elements.Norm(i, j, k)({alpha, Comp}).out;
                                double localIncrement = localDelta*(1.0 - localNorm);
                                if(localIncrement < -Precision)
                                {
                                    Elements.PhaseDot(  i,   j,   k)({alpha, Comp}).out -=
                                                                    localIncrement;
                                }
                            }
                            if (Elements.Norm(i, j, k)({alpha, Comp}).out >=
                                Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).in)
                            {
                                double localNorm = Elements.Norm(i+di, j+dj, k+dk)({alpha, Comp}).in;
                                double localIncrement = localDelta*(1.0 - localNorm);
                                if(localIncrement < -Precision)
                                {
                                    Elements.PhaseDot(  i,   j,   k)({alpha, Comp}).out -=
                                                                    localIncrement;
                                }
                            }
                        }
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    void Diffsion::ApplyIncrements(PhaseField& Phase,
                                                Composition& Elements, double dt)
    {
        double ompTotalMass = 0.0;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Elements.Total,0,reduction(+:ompTotalMass))
        {
            if(Phase.Interface(i,j,k))
            {
                for(int alpha = 0; alpha != Nphases; ++alpha)
                if(!Stoichiometric[alpha])
                {
                    double dPhase = (Elements.PhaseDot(i,j,k)({alpha, Comp}).in +
                                    Elements.PhaseDot(i,j,k)({alpha, Comp}).out)*dt;
                    if(fabs(dPhase) <= Precision) dPhase = 0.0;
                    Elements.Total(i,j,k)({Comp}) += dPhase;
                }
            }
            else
            {
                int alpha = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
                if(!Stoichiometric[alpha])
                {
                    double dPhase = (Elements.PhaseDot(i,j,k)({alpha, Comp}).in +
                                    Elements.PhaseDot(i,j,k)({alpha, Comp}).out)*dt;
                    if(fabs(dPhase) <= Precision) dPhase = 0.0;
                    Elements.Total(i,j,k)({Comp}) += dPhase;
                }
            }
            ompTotalMass += Elements.Total(i,j,k)({Comp});
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        TotalMass = ompTotalMass;
    }

    void Diffsion::Solve(PhaseField& Phase,
                        Composition& Elements, Temperature& Tx,
                        BoundaryConditions& BC, double dt,
                        bool InternalDiffusivities)
    {
        Elements.SetBoundaryConditions(BC);

        CalculatePhaseConcentrations(Phase, Tx, Elements);

        if(InternalDiffusivities)
        {
            SetDiffusionCoefficients(Phase, Tx);
        }
        CalculateDiffusionIncrements(Phase, Elements, Tx);
        //CalculateAntitrappingIncrements(Phase, Elements);
        CalculateLimits(Phase, Elements, dt);

        if(LimitingNeeded)
        {
            Elements.SetLimitsBoundaryConditions(BC);
            LimitDiffusionIncrements(Phase, Elements, Tx);
            //LimitAntitrappingIncrements(Phase, Elements, Tx);
        }
        ApplyIncrements(Phase, Elements, dt);

        Elements.SetBoundaryConditions(BC);

        CalculatePhaseConcentrations(Phase, Tx, Elements);
        if(ThereAreStoichiometricPhases)
        {
            RestoreStoichiometric(Phase, Tx, Elements);
            Elements.SetBoundaryConditions(BC);
        }
    }

    void Diffsion::PrintPointStatistics(int x, int y, int z)
    {
        cout << "Diffusion Coefficients:" << endl;
        for (int n = 0; n < Nphases; n++)
        {
            cout << "D[" << n << "] = " << DC(x, y, z)({n}) << endl;
        }
        cout << endl;
    }

    double Diffsion::ReportMaximumTimestep(PhaseField& Phase)
    {
        double maxDC = 0.0;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, DC,0,reduction(max:maxDC))
        {
            for (int n = 0; n < Nphases; ++n)
            {
                maxDC = std::max(maxDC, DC(i, j, k)({n}));
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        return 0.25*dx*dx/maxDC;
    }


}