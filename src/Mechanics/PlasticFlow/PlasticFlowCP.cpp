

#include "Mechanics/Storages/PlasticProperties.h"
#include "Info.h"
#include "Mechanics/PlasticFlow/PlasticFlowMethods.h"
#include "Mechanics/PlasticFlow/PlasticFlowCP.h"
#include "PhaseField.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "Orientations.h"
#include "BoundaryConditions.h"
#include "DrivingForce.h"
#include "Velocities.h"
#include "Chemistry/ChemicalProperties.h"
#include "Compositions.h"

namespace opensim
{
using namespace std;

PlasticFlowCP::PlasticFlowCP(Settings& locSettings)
{
    Initialize(locSettings);
}

void PlasticFlowCP::Initialize(Settings& locSettings)
{
    thisclassname = "PlasticFlowCP";
    //DefaultInputFileName = ProjectInputDir + "PlasticFlowInput.opi";
    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    Nphases = locSettings.Nphases;
    nSlipSystems = 12;                                        ///preliminary value!!!

    dt = locSettings.dt;
    HardeningSwitch = true;
    FirstStepHardening = true;
    PlasticitySwitch = false;
    TotalCRSSHardening = 0.0;
    int maxNslip = 25;
    phaseExponent.Allocate({Nphases});
    phaseSlipRate.Allocate({Nphases});
    PlasticitySwitchPhase.Allocate({Nphases});
    Nslip.Allocate({Nphases});
    BurgersLength.Allocate({Nphases});
    Hardening.Allocate({Nphases});
    phaseHardExponent.Allocate({Nphases});
    PhaseSaturationCRSS.Allocate({Nphases, maxNslip});
    InitialHardeningModulus.Allocate({Nphases, maxNslip});
    LatentHardeningParameter.Allocate({Nphases, maxNslip});
    PhaseCriticalResolvedShearStress.Allocate({Nphases, maxNslip});

    PhaseGSNormal.Allocate(Nphases);
    PhaseGSDirection.Allocate(Nphases);
    PhaseGSLine.Allocate(Nphases);
    CRSS.Allocate(Nx, Ny, Nz, 1);
    for(int alpha = 0; alpha < Nphases; alpha++)
    {
        Nslip({alpha}) = 0;
        phaseExponent({alpha}) = 0;
        phaseSlipRate({alpha}) = 0;
        phaseHardExponent({alpha}) = 0;
        Hardening({alpha}) = false;
        PhaseGSNormal[alpha].Allocate({12,3});
        PhaseGSDirection[alpha].Allocate({12,3});
        PhaseGSLine[alpha].Allocate({12,3});
        for(int slipSys = 0; slipSys < maxNslip; slipSys++)
        {
            PhaseCriticalResolvedShearStress({alpha, slipSys}) = 0.;
            PhaseSaturationCRSS({alpha, slipSys}) = 0.;
            InitialHardeningModulus({alpha, slipSys}) = 0.;
            LatentHardeningParameter({alpha, slipSys}) = 0.;
        }
    }
    initialized = true;
    Info::WriteLine();
    Info::WriteStandard(thisclassname, "Initialized");
}

void PlasticFlowCP::ReadInput(string InputFileName)
{
    thisclassname = "PlasticFlowCP";
    Info::WriteLineInsert("Plastic flow");
    Info::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened",
                                                    thisclassname, "ReadInput");
        exit(1);
    };
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
    /*stringstream namePlasticityStartStep;
    namePlasticityStartStep << "PlasticityStartStep";
    PlasticityStartStep = UserInterface::ReadParameterI(inp, moduleLocation, namePlasticityStartStep.str());*/

    for(int alpha = 0; alpha < Nphases; alpha++)
    {
        stringstream namePlasticitySwitch;
        namePlasticitySwitch << "PlasticityPhase" << alpha;
        PlasticitySwitchPhase({ alpha }) = UserInterface::ReadParameterB(inp, moduleLocation, namePlasticitySwitch.str());

        if(PlasticitySwitchPhase({ alpha }) )
        {
            PlasticitySwitch = true;
            stringstream nameCrystal;
            nameCrystal<<"Crystal"<<alpha;
            string tmp = UserInterface::ReadParameterF(inp, moduleLocation, nameCrystal.str());
            if (tmp == "FCC")
            {
                dVector3 GSNormal;
                dVector3 GSDirection;
                dVector3 GSLine;
                double LatticeStruct[9*12] =
                { //  {  n1,  n2,  n3,  d1,  d2,  d3,  l1,  l2,  l3 }
                        -1.0, 1.0, 1.0, 0.0, 1.0,-1.0,-2.0,-1.0,-1.0,
                        -1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 2.0,-1.0,
                        -1.0, 1.0, 1.0, 1.0, 1.0, 0.0,-1.0, 1.0,-2.0,

                        1.0, 1.0, 1.0, 0.0, 1.0,-1.0,-2.0, 1.0, 1.0,
                        1.0, 1.0, 1.0, 1.0, 0.0,-1.0,-1.0, 2.0,-1.0,
                        1.0, 1.0, 1.0, 1.0,-1.0, 0.0, 1.0, 1.0,-2.0,

                        1.0, 1.0,-1.0, 0.0, 1.0, 1.0, 2.0,-1.0, 1.0,
                        1.0, 1.0,-1.0, 1.0, 0.0, 1.0, 1.0,-2.0,-1.0,
                        1.0, 1.0,-1.0, 1.0,-1.0, 0.0,-1.0,-1.0,-2.0,

                        1.0,-1.0, 1.0, 0.0, 1.0, 1.0,-2.0,-1.0, 1.0,
                        1.0,-1.0, 1.0, 1.0, 0.0,-1.0, 1.0, 2.0, 1.0,
                        1.0,-1.0, 1.0, 1.0, 1.0, 0.0,-1.0, 1.0, 2.0
                };
                for(int ss = 0; ss < 12; ss ++)
                {
                    for(int j = 0; j < 3; ++j)
                    {
                        GSNormal[j]    = LatticeStruct[ss*9+j];
                        GSDirection[j] = LatticeStruct[ss*9+j+3];
                        GSLine[j]  = LatticeStruct[ss*9+j+6];
                    }
                    GSNormal.normalize();
                    GSDirection.normalize();
                    GSLine.normalize();

                    for(int j = 0; j < 3; ++j)
                    {
                        PhaseGSDirection[alpha]({ss,j}) = GSDirection[j];
                        PhaseGSNormal[alpha]({ss,j}) = GSNormal[j];
                        PhaseGSLine[alpha]({ss,j}) = GSLine[j];
                    }
                }
            }
            else if (tmp == "BCC")
            {
                dVector3 GSNormal;
                dVector3 GSDirection;
                dVector3 GSLine;
                double LatticeStruct[9*12] =
                { //  {  n1,  n2,  n3,  d1,  d2,  d3,  l1,  l2,  l3 }
                        0.0, 1.0,-1.0, 1.0, 1.0, 1.0, 2.0,-1.0,-1.0,
                        1.0, 0.0,-1.0, 1.0, 1.0, 1.0, 1.0,-2.0, 1.0,
                        1.0,-1.0, 0.0, 1.0, 1.0, 1.0,-1.0,-1.0, 2.0,

                        0.0, 1.0,-1.0,-1.0, 1.0, 1.0,-2.0, 1.0, 1.0,
                        1.0, 0.0, 1.0,-1.0, 1.0, 1.0,-1.0, 2.0, 1.0,
                        1.0, 1.0, 0.0,-1.0,-1.0, 1.0, 1.0,-1.0, 2.0,

                        0.0, 1.0, 1.0, 1.0, 1.0,-1.0,-2.0, 1.0,-1.0,
                        1.0, 0.0, 1.0, 1.0, 1.0, 1.0,-1.0, 2.0, 1.0,
                        1.0,-1.0, 0.0, 1.0, 1.0,-1.0, 1.0, 1.0, 2.0,

                        0.0, 1.0, 1.0, 1.0,-1.0, 1.0, 2.0, 1.0,-1.0,
                        1.0, 0.0,-1.0, 1.0,-1.0, 1.0, 1.0, 2.0, 1.0,
                        1.0, 1.0, 0.0, 1.0,-1.0, 1.0, 1.0,-1.0,-2.0
                };
                for(int ss = 0; ss < 12; ss ++)
                {
                    for(int j = 0; j < 3; ++j)
                    {
                        GSNormal[j]    = LatticeStruct[ss*9+j];
                        GSDirection[j] = LatticeStruct[ss*9+j+3];
                        GSLine[j]  = LatticeStruct[ss*9+j+6];
                    }
                    GSNormal.normalize();
                    GSDirection.normalize();
                    GSLine.normalize();

                    for(int j = 0; j < 3; ++j)
                    {
                        PhaseGSDirection[alpha]({ss,j}) = GSDirection[j];
                        PhaseGSNormal[alpha]({ss,j}) = GSNormal[j];
                        PhaseGSLine[alpha]({ss,j}) = GSLine[j];
                    }
                }
            }
            else
            {
                Info::WriteExit("Crystal structure invalid",thisclassname, "ReadInput()");
                exit(1);
            }

            stringstream nameNslip;
            nameNslip << "Nslip" << alpha;
            Nslip({alpha}) = 12;

            stringstream namePhaseExponent;
            namePhaseExponent << "n" << alpha;
            phaseExponent({alpha}) = UserInterface::ReadParameterD(inp, moduleLocation, namePhaseExponent.str());

            stringstream namePhaseSlipRate;
            namePhaseSlipRate << "Sliprate" << alpha;
            phaseSlipRate({alpha}) = UserInterface::ReadParameterD(inp, moduleLocation, namePhaseSlipRate.str());

            stringstream nameHardening;
            nameHardening << "Hardening" << alpha;
            Hardening({alpha}) = UserInterface::ReadParameterB(inp, moduleLocation, nameHardening.str());

            stringstream namePhaseHardEponent;
            namePhaseHardEponent << "HardEponent" << alpha;
            phaseHardExponent({alpha}) = UserInterface::ReadParameterD(inp, moduleLocation, namePhaseHardEponent.str());

            stringstream nameBurgers;
            nameBurgers << "Burgers" << alpha;
            BurgersLength({alpha}) =
                    UserInterface::ReadParameterD(inp, moduleLocation, nameBurgers.str());

            for(int slipSys = 0; slipSys < Nslip({alpha}); slipSys++)
            {
                stringstream nameCRSS;
                nameCRSS<<"CRSS"<<alpha<<"_"<< slipSys;
                PhaseCriticalResolvedShearStress({alpha, slipSys}) =
                        UserInterface::ReadParameterD(inp, moduleLocation, nameCRSS.str());
            }
            for(int slipSys = 0; slipSys < Nslip({alpha}); slipSys++)
            {
                stringstream namesCRSS;
                namesCRSS<<"sCRSS"<<alpha<<"_"<< slipSys;
                PhaseSaturationCRSS({alpha, slipSys}) =
                        UserInterface::ReadParameterD(inp, moduleLocation, namesCRSS.str());
            }
            for(int slipSys = 0; slipSys < Nslip({alpha}); slipSys++)
            {
                stringstream nameHard;
                nameHard<<"IniHard"<<alpha<<"_"<< slipSys;
                InitialHardeningModulus({alpha, slipSys}) =
                        UserInterface::ReadParameterD(inp, moduleLocation, nameHard.str());
            }
            for(int slipSys = 0; slipSys < Nslip({alpha}); slipSys++)
            {
                stringstream nameLatent;
                nameLatent<<"Latent"<<alpha<<"_"<< slipSys;
                LatentHardeningParameter({alpha, slipSys}) =
                        UserInterface::ReadParameterD(inp, moduleLocation, nameLatent.str());
            }
            Info::WriteLine();
        }
    }

	stringstream nameCutoff;
	nameCutoff << "ShearRateCutoff";
	allowedShearRate = UserInterface::ReadParameterD(inp, moduleLocation, nameCutoff.str());

	Info::WriteLine();
}

void PlasticFlowCP::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(CRSS);
    BC.SetY(CRSS);
    BC.SetZ(CRSS);
}

void PlasticFlowCP::Remesh(int newNx, int newNy, int newNz,
                                                         BoundaryConditions& BC)
{
    SetBoundaryConditions(BC);

    CRSS.Remesh(newNx, newNy, newNz);

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    SetBoundaryConditions(BC);
    Info::WriteStandard(thisclassname, "Remeshed");
}

void PlasticFlowCP::Advect(PhaseField& Phase, Velocities& Vel, BoundaryConditions& BC,
                    double dt, int scheme)
{
    if(CRSSdot.IsNotAllocated())
    {
        CRSSdot.Allocate(Nx, Ny, Nz, 1);
    }
    if(not CRSSdot.IsSize(Nx, Ny, Nz))
    {
        CRSSdot.Reallocate(Nx, Ny, Nz);
    }
    SetBoundaryConditions(BC);

    switch(scheme)
    {
        case Upwind:
        {
            const double dx = Vel.dx;
            const double dx2 = 0.5/dx;

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,CRSSdot,0,)
            {
                double ux = Vel.Average(i,j,k)[0];
                double uy = Vel.Average(i,j,k)[1];
                double uz = Vel.Average(i,j,k)[2];

                double uxm = Vel.Average(i-1,j,k)[0];
                double uxp = Vel.Average(i+1,j,k)[0];

                double uym = Vel.Average(i,j-1,k)[1];
                double uyp = Vel.Average(i,j+1,k)[1];

                double uzm = Vel.Average(i,j,k-1)[2];
                double uzp = Vel.Average(i,j,k+1)[2];


                std::vector<int> tempPFindex = Phase.ReturnVicinityPhaseFields(i,j,k);

                for(auto it = tempPFindex.begin(); it != tempPFindex.end(); ++it)
                {
                    int thPhaseIndex = Phase.FieldsStatistics[*it].Phase;
                    if (Hardening[thPhaseIndex])
                    {
                        for (int ss = 0; ss < 12; ss++)
                        {
                            double ijk = 0;
                            if(Phase.PhaseFieldPresent(i,j,k,*it))
                            {
                                    ijk = CRSS(i,j,k).get(*it, ss);
                            }
                            else
                            {
                                ijk = 0.0;
                            }

                            //i+1 j k
                            double i1jk = 0;
                            if(Phase.PhaseFieldPresent(i+1,j,k,*it))
                            {
                                    i1jk = CRSS(i+1,j,k).get(*it, ss);
                            }
                            else
                            {
                                i1jk = 0.0;
                            }

                            //i-1 j k
                            double im1jk = 0;
                            if(Phase.PhaseFieldPresent(i-1,j,k,*it))
                            {
                                    im1jk = CRSS(i-1,j,k).get(*it, ss);
                            }
                            else
                            {
                                im1jk = 0.0;
                            }

                            //i j+1 k
                            double ij1k = 0;
                            if(Phase.PhaseFieldPresent(i,j+1,k,*it))
                            {
                                    ij1k = CRSS(i,j+1,k).get(*it, ss);
                            }
                            else
                            {
                                ij1k = 0.0;
                            }

                            //i j-1 k
                            double ijm1k = 0;
                            if(Phase.PhaseFieldPresent(i,j-1,k,*it))
                            {
                                    ijm1k = CRSS(i,j-1,k).get(*it, ss);
                            }
                            else
                            {
                                ijm1k = 0.0;
                            }

                            //i j k+1
                            double ijk1 = 0;
                            if(Phase.PhaseFieldPresent(i,j,k+1,*it))
                            {
                                    ijk1 = CRSS(i,j,k+1).get(*it, ss);
                            }
                            else
                            {
                                ijk1 = 0.0;
                            }

                            //i j k-1
                            double ijkm1 = 0;
                            if(Phase.PhaseFieldPresent(i,j,k-1,*it))
                            {
                                    ijkm1 = CRSS(i,j,k-1).get(*it, ss);
                            }
                            else
                            {
                                ijkm1 = 0.0;
                            }

                            CRSSdot(i,j,k).set(*it, ss,
                                (im1jk*(fabs(uxm) + uxm) +
                                 ijm1k*(fabs(uym) + uym) +
                                 ijkm1*(fabs(uzm) + uzm) +
                                 i1jk*(fabs(uxp) - uxp) +
                                 ij1k*(fabs(uyp) - uyp) +
                                 ijk1*(fabs(uzp) - uzp))*dx2 -
                                 ijk*(fabs(ux) + fabs(uy) + fabs(uz))*(1.0/dx));
                        }
                    }
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            break;
        }
        case LaxWendroff:
        {
            Info::WriteExit("LaxWendroff advection scheme not supported for CRSS",
                            thisclassname, "AdvectPlasticStrain(PFP, Vel, BC, dt)");
            exit(13);
            break;
        }
        default:
        {
            Info::WriteExit("Wrong/None advection scheme given in the input file",
                             thisclassname, "AdvectPlasticStrain(PFP, Vel, BC, dt)");
            exit(13);
        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,CRSSdot,0,)
    {
        for (auto it = Phase.Fields(i,j,k).cbegin();
                  it < Phase.Fields(i,j,k).cend(); ++it)
        {
            int phaseIndex = it->index;
            int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;

            for(int slipSys = 0; slipSys < Nslip({thPhaseIndex}); slipSys++)
            {
                double CRSSOld = CRSS(i,j,k).get(phaseIndex, slipSys);
                double temp = CRSSOld + CRSSdot(i, j, k).get(phaseIndex, slipSys) * dt;
                double satCRSS = PhaseSaturationCRSS({thPhaseIndex, slipSys});
                if (temp > satCRSS)
                    temp = satCRSS;
                CRSS(i, j, k).set(phaseIndex, slipSys, temp);
                CRSSdot(i, j, k).set(phaseIndex, slipSys, 0.0);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    SetBoundaryConditions(BC);
}

double PlasticFlowCP::Solve(vector<OPObject*> Objects, bool verbose)
{
//    PhaseField& Phase           = *(static_cast<PhaseField*>(findOPObject(Objects, "PhaseField", thisclassname, "Solve", true)));
//    BoundaryConditions& BC      = *(static_cast<BoundaryConditions*>(findOPObject(Objects, "BoundaryConditions", thisclassname, "Solve", true)));
//    ElasticProperties& EP       = *(static_cast<ElasticProperties*>(findOPObject(Objects, "ElasticProperties", thisclassname, "Solve", true)));
//    Orientations& OR            = *(static_cast<Orientations*>(findOPObject(Objects, "Orientations", thisclassname, "Solve", true)));
//    PlasticFlProperties& PFP  = *(static_cast<PlasticFlProperties*>(findOPObject(Objects, "PlasticFlProperties", thisclassname, "Solve", true)));
	return PlasticFlowMethods::CalculatePlasticStrainCP(Objects, verbose);
}

void PlasticFlowCP::SetInitialHardening(const PhaseField& Phase, BoundaryConditions& BC)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,CRSS,0,)
    {
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin(); alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                int phaseIndex = alpha->index;
                int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
//                if(Hardening[thPhaseIndex])
//                {
                    for(int slipSys = 0; slipSys < 12; slipSys++)
                    {
                        double iniCRSS = PhaseCriticalResolvedShearStress({thPhaseIndex, slipSys});
                        CRSS(i,j,k).set(phaseIndex, slipSys, iniCRSS);
                    }
//                }
            } // end alpha
        }
        else
        {
            int phaseIndex = Phase.Fields(i,j,k).front().index;
            int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
//            if(Hardening[thPhaseIndex])
//            {
                for(int slipSys = 0; slipSys < 12; slipSys++)
                {
                    double iniCRSS = PhaseCriticalResolvedShearStress({thPhaseIndex, slipSys});
                    CRSS(i,j,k).set(phaseIndex, slipSys, iniCRSS);
                }
//            }
        } // end interface
    } // end loop space
    OMP_PARALLEL_STORAGE_LOOP_END
    SetBoundaryConditions(BC);
}

void PlasticFlowCP::SetInitialHardening(const PhaseField& Phase, BoundaryConditions& BC, int targetPhaseIndex)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,CRSS,0,)
    {
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin(); alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                int phaseIndex = alpha->index;
                int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
                if(targetPhaseIndex == phaseIndex)
                {
                    for(int slipSys = 0; slipSys < 12; slipSys++)
                    {
                        double iniCRSS = PhaseCriticalResolvedShearStress({thPhaseIndex, slipSys});
                        CRSS(i,j,k).set(phaseIndex, slipSys, iniCRSS);
                    }
                }
            } // end alpha
        }
        else
        {
            int phaseIndex = Phase.Fields(i,j,k).front().index;
            int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
            if(targetPhaseIndex == phaseIndex)
            {
                for(int slipSys = 0; slipSys < 12; slipSys++)
                {
                    double iniCRSS = PhaseCriticalResolvedShearStress({thPhaseIndex, slipSys});
                    CRSS(i,j,k).set(phaseIndex, slipSys, iniCRSS);
                }
            }
        } // end interface
    } // end loop space
    OMP_PARALLEL_STORAGE_LOOP_END
    SetBoundaryConditions(BC);
}

void PlasticFlowCP::SetInitialHardening(const PhaseField& Phase, const Composition& Cx, BoundaryConditions& BC)
{
    double m = 3.69713e+10;
    double b = 7.48763e+08;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,CRSS,0,)
    {
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin(); alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                int phaseIndex = alpha->index;
                int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
                for(int slipSys = 0; slipSys < 12; slipSys++)
                {
                    if (thPhaseIndex == 1)
                    {
//                        double b = PhaseCriticalResolvedShearStress({thPhaseIndex, slipSys});
                        double iniCRSS = m * Cx.Phase(i,j,k)({thPhaseIndex, 0}) + b;
                        CRSS(i,j,k).set(phaseIndex, slipSys, iniCRSS);
                    }
                    else
                    {
                        double iniCRSS = PhaseCriticalResolvedShearStress({thPhaseIndex, slipSys});
                        CRSS(i,j,k).set(phaseIndex, slipSys, iniCRSS);
                    }
                }
            } // end alpha
        }
        else
        {
            int phaseIndex = Phase.Fields(i,j,k).front().index;
            int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
            for(int slipSys = 0; slipSys < 12; slipSys++)
            {
                if (thPhaseIndex == 1)
                {
//                    double b = PhaseCriticalResolvedShearStress({thPhaseIndex, slipSys});
                    double iniCRSS = m * Cx.Phase(i,j,k)({thPhaseIndex, 0}) + b;
                    CRSS(i,j,k).set(phaseIndex, slipSys, iniCRSS);
                }
                else
                {
                    double iniCRSS = PhaseCriticalResolvedShearStress({thPhaseIndex, slipSys});
                    CRSS(i,j,k).set(phaseIndex, slipSys, iniCRSS);
                }
            }
        } // end interface
    } // end loop space
    OMP_PARALLEL_STORAGE_LOOP_END
    SetBoundaryConditions(BC);
    for(int slipSys = 0; slipSys < 12; slipSys++)
    {
        PhaseSaturationCRSS({1, slipSys}) = (m * 0.00497673 + b) * 1.05;
    }
}

dVector6 PlasticFlowCP::GetSchmidMatrixSym(int i, int j, int k, int phaseIndex, int slipSys, PhaseField& Phase, ElasticProperties& EP, Orientations& OR)
{
    int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
	int vIndex = Phase.FieldsStatistics[thPhaseIndex].Variant;
    dVector6 tmp;
    dVector3 GSDirection;
    dVector3 GSNormal;
    for (int i=0; i<3; i++)
    {
        GSDirection[i] = PhaseGSDirection[thPhaseIndex]({slipSys,i});
        GSNormal[i] = PhaseGSNormal[thPhaseIndex]({slipSys,i});
    }

	if (EP.Variants.set)
	{
		GSDirection.rotate(EP.Variants(phaseIndex, vIndex));
		GSNormal.rotate(EP.Variants(phaseIndex, vIndex));
	}

    dMatrix3x3 tempRotMatrix = OR.Quaternions(i,j,k).RotationMatrix * Phase.FieldsStatistics[phaseIndex].Orientation.RotationMatrix;
    GSDirection = GSDirection.rotated(tempRotMatrix);
    GSNormal = GSNormal.rotated(tempRotMatrix);
    GSDirection.normalize();
    GSNormal.normalize();
    tmp[0] = GSDirection[0]*GSNormal[0];
    tmp[1] = GSDirection[1]*GSNormal[1];
    tmp[2] = GSDirection[2]*GSNormal[2];
    tmp[3] = 0.5*(GSDirection[1]*GSNormal[2] +
                                     GSDirection[2]*GSNormal[1]);
    tmp[4] = 0.5*(GSDirection[0]*GSNormal[2] +
                                     GSDirection[2]*GSNormal[0]);
    tmp[5] = 0.5*(GSDirection[0]*GSNormal[1] +
                                     GSDirection[1]*GSNormal[0]);
    return tmp;
}

dMatrix3x3 PlasticFlowCP::GetSchmidMatrix(int i, int j, int k, int phaseIndex, int slipSys, PhaseField& Phase, ElasticProperties& EP, Orientations& OR)
{
    int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
	int vIndex = Phase.FieldsStatistics[thPhaseIndex].Variant;
    dMatrix3x3 tmp;
    dVector3 GSDirection;
    dVector3 GSNormal;
    for (int i=0; i<3; i++)
    {
        GSDirection[i] = PhaseGSDirection[thPhaseIndex]({slipSys,i});
        GSNormal[i] = PhaseGSNormal[thPhaseIndex]({slipSys,i});
    }

	if (EP.Variants.set)
	{
		GSDirection.rotate(EP.Variants(phaseIndex, vIndex));
		GSNormal.rotate(EP.Variants(phaseIndex, vIndex));
	}

    dMatrix3x3 tempRotMatrix = OR.Quaternions(i,j,k).RotationMatrix * Phase.FieldsStatistics[phaseIndex].Orientation.RotationMatrix;
    GSDirection = GSDirection.rotated(tempRotMatrix);
    GSNormal = GSNormal.rotated(tempRotMatrix);
    GSDirection.normalize();
    GSNormal.normalize();
    tmp(0,0) = GSDirection[0]*GSNormal[0];
    tmp(1,1) = GSDirection[1]*GSNormal[1];
    tmp(2,2) = GSDirection[2]*GSNormal[2];

    tmp(0,1) = GSDirection[0]*GSNormal[1];
    tmp(0,2) = GSDirection[0]*GSNormal[2];
    tmp(1,2) = GSDirection[1]*GSNormal[2];

    tmp(1,0) = GSDirection[1]*GSNormal[0];
    tmp(2,0) = GSDirection[2]*GSNormal[0];
    tmp(2,1) = GSDirection[2]*GSNormal[1];

//    tmp.rotated(tempRotMatrix);

    return tmp;
}

void PlasticFlowCP::Write(int tStep, bool legacy_format)
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"CRSS_", tStep, ".dat");

    ofstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be created", thisclassname, "Write()");
        exit(1);
    };

    if(not legacy_format)
    {
        int tmp = Nx;
        out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
        tmp = Ny;
        out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
        tmp = Nz;
        out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
    }
    STORAGE_LOOP_BEGIN(i,j,k,CRSS,0)
    {
        int tmp = CRSS(i,j,k).size();
        out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
        for(auto n = CRSS(i,j,k).cbegin();
                 n < CRSS(i,j,k).cend(); ++n)
        {
            int idx = n->index;
            out.write(reinterpret_cast<char*>(&idx), sizeof(int));
            for (int ss = 0; ss < 12; ss++)
            {
                double tempCRSS = n->value[ss];
                out.write(reinterpret_cast<char*>(&tempCRSS), sizeof(double));
            }
        }
    }
    STORAGE_LOOP_END
    out.close();
}

void PlasticFlowCP::Read(const BoundaryConditions& BC, int tStep, bool legacy_format)
{
    string FileName =
        UserInterface::MakeFileName(RawDataDir,"CRSS_", tStep, ".dat");

    Read(FileName, legacy_format);
    SetBoundaryConditions(BC);
}

void PlasticFlowCP::Read(string FileName, bool legacy_format)
{
    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit(FileName + " could not be opened",
                thisclassname, "Read()");
        exit(1);
    };

    if(not legacy_format)
    {
        int locNx = Nx;
        int locNy = Ny;
        int locNz = Nz;
        inp.read(reinterpret_cast<char*>(&locNx), sizeof(int));
        inp.read(reinterpret_cast<char*>(&locNy), sizeof(int));
        inp.read(reinterpret_cast<char*>(&locNz), sizeof(int));
        if(locNx != Nx or locNy != Ny or locNz != Nz)
        {
            stringstream message;
            message << "Inconsistent system dimensions!\n"
                    << "Input data dimensions: (" << locNx
                    << ", " << locNy << ", " << locNz << ") grid points.\n"
                    << "Required data dimensions: (" << Nx
                    << ", " << Ny << ", " << Nz << ") grid points.\n";
            Info::WriteExit(message.str(), thisclassname, "Read()");
            exit(1);
        }
    }
    STORAGE_LOOP_BEGIN(i,j,k,CRSS,0)
    {
        int   num = 0;
        int   idx = 0;
        double val = 0.0;
        CRSS(i,j,k).clear();
        inp.read(reinterpret_cast<char*>(&num), sizeof(int));                   // Fields(i,j,k).size()
        for(int n = 0; n < num; n++)
        {
            inp.read(reinterpret_cast<char*>(&idx), sizeof(int));               // Fields(i,j,k)->index
            for (int ss = 0; ss < 12; ss++)
            {
                inp.read(reinterpret_cast<char*>(&val), sizeof(double));            // Fields(i,j,k)->value
                CRSS(i,j,k).set(idx, ss, val);
            }
        }
    }
    STORAGE_LOOP_END
    inp.close();

    Info::WriteStandard(thisclassname, "Binary input loaded");
}


void PlasticFlowCP::WriteCRSSVTK(PhaseField& Phase, int tStep)
{
    stringstream outbufer;

    outbufer << "# vtk DataFile Version 3.0\n";
    outbufer << "CRSS\n";
    outbufer << "ASCII\n";
    outbufer << "DATASET STRUCTURED_GRID\n";
    outbufer << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    outbufer << "POINTS " <<  Nx*Ny*Nz << " double\n";

    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        outbufer << i << " " << j << " " << k << "\n";
    }
    outbufer << "\n";
    outbufer << "POINT_DATA " << Nx*Ny*Nz << "\n";
    for (int ss = 0; ss < 12; ss++)
    {
        outbufer << "SCALARS CRSS_" << ss << " double 1" << " \n";
        outbufer << "LOOKUP_TABLE default" << " \n";
        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            double sum = 0.;
            if(Phase.Interface(i,j,k))
            for (auto alpha = Phase.Fields(i,j,k).cbegin();
                      alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                int phaseIndex = alpha->index;
                int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
                if (alpha->value != 0)
                {
                    if (Hardening({thPhaseIndex}))
                    {
                        sum += CRSS(i,j,k).get(alpha->index, ss) * alpha->value;
                    }
                    else
                    {
                        sum += PhaseCriticalResolvedShearStress({thPhaseIndex, ss}) * alpha->value;
                    }
                }
            }
            else
            {
                int phaseIndex = Phase.Fields(i,j,k).front().index;
                int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
                if (Hardening({thPhaseIndex}))
                {
                    sum = CRSS(i,j,k).get(phaseIndex, ss);
                }
                else
                {
                    sum = PhaseCriticalResolvedShearStress({thPhaseIndex, ss});
                }
                sum = CRSS(i,j,k).get(Phase.Fields(i,j,k).front().index, ss);
            }
        outbufer << sum << "\n";
        }
    }

    string FileName = UserInterface::MakeFileName(VTKDir,"CRSS_", tStep, ".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbufer.rdbuf();
    vtk_file.close();
}
}