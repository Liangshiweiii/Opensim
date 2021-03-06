

#include <Mechanics/Storages/PlasticProperties.h>
#include "Info.h"
#include "Mechanics/PlasticFlow/PlasticFlowMethods.h"
#include "Mechanics/PlasticFlow/PlasticFlow.h"
#include "Mechanics/PlasticFlow/PlasticFlowCP.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "Mechanics/Storages/DamageProperties.h"
#include "Orientations.h"
#include "BoundaryConditions.h"
#include "Initializations.h"
#include "Compositions.h"
#include "Temperatures.h"

namespace opensim
{
using namespace std;

double PlasticFlowMethods::CalculatePlasticStrainCP(std::vector<OPObject*> Objects, bool verbose)
{
    string fctName = "CalculatePlasticStrain";
    string className = "PlasticFlowMethods";

    PhaseField& Phase           = *(static_cast<PhaseField*>(OPObject::findOPObject(Objects, "PhaseField", className, fctName, true)));
    Orientations& OR            = *(static_cast<Orientations*>(OPObject::findOPObject(Objects, "Orientations", className, fctName, true)));
    PlasticProperties& PFP    = *(static_cast<PlasticProperties*>(OPObject::findOPObject(Objects, "PlasticProperties", className, fctName, true)));
    ElasticProperties& EP       = *(static_cast<ElasticProperties*>(OPObject::findOPObject(Objects, "ElasticProperties", className, fctName, true)));
    PlasticFlowCP& PF           = *(static_cast<PlasticFlowCP*>(OPObject::findOPObject(Objects, "PlasticFlowCP", className, fctName, true)));
    BoundaryConditions& BC      = *(static_cast<BoundaryConditions*>(OPObject::findOPObject(Objects, "BoundaryConditions", className, fctName, true)));
    Composition* Cx             = (static_cast<Composition*>(OPObject::findOPObject(Objects, "Composition", className, fctName, false, false)));
    DamageProperties* DP        = (static_cast<DamageProperties*>(OPObject::findOPObject(Objects, "DamageProperties", className, fctName, false, false)));

    double dt = PF.dt;
    double maxStrainNorm = 0.0;
    double allowedShearRate = PF.allowedShearRate;
    double MAXSROvershoot = 0.0;
    int sRLimitWarning = 0;
    bool cxCRSS = false;
    if (Cx != nullptr)
    {
        cxCRSS = true;
    }
    bool damage = false;
    if (DP != nullptr)
    {
        damage = true;
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields, 0,)
    {
        vStrain plasticStrainIncrement;
        plasticStrainIncrement.set_to_zero();

        for (auto it = Phase.Fields(i,j,k).cbegin();
                  it != Phase.Fields(i,j,k).cend(); ++it)
        {
            if (it->value != 0)
            {
                int phaseIndex = it->index;
                int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;

                dVector<12> shearRate;
                shearRate.set_to_zero();
                if(PF.PlasticitySwitchPhase({thPhaseIndex}));
                for(int slipSys = 0; slipSys < PF.Nslip({thPhaseIndex}); slipSys++)
                {
                    dMatrix3x3 locSchmidMatrix =
                                    PF.GetSchmidMatrix(i,j,k, phaseIndex, slipSys, Phase, EP, OR);

                    dMatrix3x3 locStress= EP.Stresses(i,j,k).tensor();

                    double resolvedShear = locStress.double_contract(locSchmidMatrix);
                    double currCRSS = 0;
                    if (PF.Hardening({thPhaseIndex}) and PF.HardeningSwitch)
                    {
                        currCRSS = PF.CRSS(i,j,k).get(phaseIndex, slipSys);
                    }
                    else
                    {
                        currCRSS = PF.PhaseCriticalResolvedShearStress({thPhaseIndex, slipSys});
                    }
                    if (damage)
                    {
                        currCRSS *= (1. - DP->EffectiveDamage(i,j,k));
                    }

                    if(currCRSS == 0.0)
                    	 cout << "currCRSS is 0 ("<< i << "," << j << "," << k << ")"   << endl;

                    double tmpShearRate = PF.phaseSlipRate({thPhaseIndex})*
                                                     pow(fabs(resolvedShear/currCRSS),
                                                  PF.phaseExponent({thPhaseIndex}));

                    if (tmpShearRate > 0.4*allowedShearRate)
                    {
                        sRLimitWarning++;
                        MAXSROvershoot = max(MAXSROvershoot, tmpShearRate / (0.4 * allowedShearRate));
                    }
                    shearRate[slipSys] = fabs(allowedShearRate * tanh(tmpShearRate / allowedShearRate));

                    maxStrainNorm = max(maxStrainNorm, fabs(shearRate[slipSys]));
                    if(resolvedShear < 0)
                        shearRate[slipSys] *= -1;
                    dMatrix3x3 tmpStrainTensor = locSchmidMatrix.getsym();
                    vStrain tmpStrain = tmpStrainTensor.VoigtStrain();
                    tmpStrain *= it->value*shearRate[slipSys]*dt;
                    plasticStrainIncrement += tmpStrain;

                }
                if (PF.Hardening({thPhaseIndex}) and PF.HardeningSwitch)
                {
                    if (!cxCRSS)
                    {
                        ApplyHardening(Phase, PF, PFP, shearRate, phaseIndex, i, j, k);
                    }
                    else
                    {
                        ApplyHardening(Phase, PF, PFP, *Cx, shearRate, phaseIndex, i, j, k);
                    }
                }
            }
        }

        vStrain currAppliedStrain;
        if (EP.LargeDeformations)
        {
			currAppliedStrain =  EP.StrainIncrements(i, j, k).Ln();//;EP.EffectiveAppliedStrain.Ln() * 1.0;
        }
        else
        {
            currAppliedStrain= EP.Strains(i,j,k);
        }
        double maxAppliedStrain = DBL_EPSILON;
        int principalDir = 0;
        for (int n = 0; n < 6; n++)
        {
             if(fabs(currAppliedStrain[n]) > maxAppliedStrain)
             {
                 maxAppliedStrain = fabs(currAppliedStrain[n]);
                 principalDir = n;
             }
        }
		double corrFactor = 1.0;
		double absCurAppliedStrain = fabs(currAppliedStrain[principalDir]);
		if (absCurAppliedStrain > 0)
		corrFactor = fabs(plasticStrainIncrement[principalDir]) / absCurAppliedStrain;
        if(corrFactor > 1.0)
        {
			 //cout << corrFactor << " = " << fabs(plasticStrainIncrement[principalDir]) << " / " << fabs(currAppliedStrain[principalDir]) << endl;
             for (int n = 0; n < 6; n++)
             {
                 plasticStrainIncrement[n] /= corrFactor;
             }
        }

        PFP.PlasticStrain(i,j,k) += plasticStrainIncrement;
//        maxStrainNorm = max(plasticStrainIncrement.norm(), maxStrainNorm);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    if (sRLimitWarning and PF.HardeningSwitch and verbose)
    {
        std::string message = "The shearrate value has been limited " +
                               std::to_string(static_cast <long int> (sRLimitWarning)) + " times\n"
                               + "          Max overshoot : " + std::to_string(static_cast <long double> (MAXSROvershoot)) + " times\n";
        Info::WriteWarning(message, "PlasticFlowCP", "PrintDiagnostics()", false);
    }
    PFP.SetBoundaryConditions(BC);
    return maxStrainNorm;///(Nx*Ny*Nz);
}

void PlasticFlowMethods::ApplyHardening(PhaseField& Phase, PlasticFlowCP& PF, PlasticProperties& PFP, dVector<12> shearRate, int phaseIndex, int i, int j, int k)
{
    int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
    double dt = PF.dt;
    int Nslip = PF.Nslip({thPhaseIndex});

//    double phaseStrain = 0.0;
    vector<double> HardeningModulus(Nslip);
    double incCRSS = 0.0;

    for(int slipSys = 0; slipSys < Nslip; slipSys++)
    {
//        phaseStrain += fabs(shearRate[slipSys])*dt;
        HardeningModulus[slipSys] = 0.0;
    }

    for(int slipSys = 0; slipSys < Nslip; slipSys++)
    {
        double iniHardModulus = PF.InitialHardeningModulus({thPhaseIndex, slipSys});
        double currCRSS = PF.CRSS(i,j,k).get(phaseIndex, slipSys);
        double satCRSS = PF.PhaseSaturationCRSS({thPhaseIndex, slipSys});
        double tauRatio = fabs(currCRSS/satCRSS);
        if(tauRatio < 1)
        {
            HardeningModulus[slipSys] = iniHardModulus*pow((1-tauRatio), PF.phaseHardExponent({thPhaseIndex}));
        }
        else
        {
            HardeningModulus[slipSys] = 0.0;
        }
    }

    for(int slipSys = 0; slipSys < Nslip; slipSys++)
    {
        double hardModulus = HardeningModulus[slipSys];
        double gammaDot = fabs(shearRate[slipSys]);
        double currCRSS = PF.CRSS(i,j,k).get(phaseIndex, slipSys);
        double satCRSS = PF.PhaseSaturationCRSS({thPhaseIndex, slipSys});
        incCRSS = hardModulus*gammaDot*dt;
        for(int jSlipSys = 0; jSlipSys< Nslip; jSlipSys++)
        {
            if(slipSys != jSlipSys)
            {
                incCRSS += PF.LatentHardeningParameter({thPhaseIndex, jSlipSys})*
                        hardModulus*fabs(shearRate[jSlipSys]) * dt;// * Phase.Fields(i,j,k)[phaseIndex];
            }
        }
        double newCRSS = currCRSS + incCRSS;
        if (newCRSS > satCRSS)
        {
            newCRSS = satCRSS;
        }
        PF.CRSS(i,j,k).set(phaseIndex, slipSys, newCRSS);
    }
}

void PlasticFlowMethods::ApplyHardening(PhaseField& Phase, PlasticFlowCP& PF, PlasticProperties& PFP, Composition& Cx, dVector<12> shearRate, int phaseIndex, int i, int j, int k)
{
    int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
    double dt = PF.dt;
    int Nslip = PF.Nslip({thPhaseIndex});
    double m = 10.0;
    double b = 12.0e8;

//    double phaseStrain = 0.0;
    vector<double> HardeningModulus(Nslip);
    double incCRSS = 0.0;

    for(int slipSys = 0; slipSys < Nslip; slipSys++)
    {
//        phaseStrain += fabs(shearRate[slipSys])*dt;
        HardeningModulus[slipSys] = 0.0;
    }

    for(int slipSys = 0; slipSys < Nslip; slipSys++)
    {
        double iniHardModulus = PF.InitialHardeningModulus({thPhaseIndex, slipSys});
        double currCRSS = PF.CRSS(i,j,k).get(phaseIndex, slipSys);
        double satCRSS = m * Cx.Phase(i,j,k)({thPhaseIndex, 0}) + b;
        double tauRatio = fabs(currCRSS/satCRSS);
        if(tauRatio < 1)
        {
            HardeningModulus[slipSys] = iniHardModulus*pow((1-tauRatio), PF.phaseHardExponent({thPhaseIndex}));
        }
        else
        {
            HardeningModulus[slipSys] = 0.0;
        }
    }

    for(int slipSys = 0; slipSys < Nslip; slipSys++)
    {
        double hardModulus = HardeningModulus[slipSys];
        double gammaDot = fabs(shearRate[slipSys]);
        double currCRSS = PF.CRSS(i,j,k).get(phaseIndex, slipSys);
        double satCRSS = PF.PhaseSaturationCRSS({thPhaseIndex, slipSys});
        incCRSS = hardModulus*gammaDot*dt;
        for(int jSlipSys = 0; jSlipSys< Nslip; jSlipSys++)
        {
            if(slipSys != jSlipSys)
            {
                incCRSS += PF.LatentHardeningParameter({thPhaseIndex, jSlipSys})*
                        hardModulus*fabs(shearRate[jSlipSys]) * dt;// * Phase.Fields(i,j,k)[phaseIndex];
            }
        }
        double newCRSS = currCRSS + incCRSS;
        if (newCRSS > satCRSS)
        {
            newCRSS = satCRSS;
        }
        PF.CRSS(i,j,k).set(phaseIndex, slipSys, newCRSS);
    }
}

void PlasticFlowMethods::UpdateCRSS(const PhaseField& Phase, PlasticFlowCP& PFCP, BoundaryConditions& BC)
{
    Tensor<double, 2> avgCRSS;
    int numberOfGrains = Phase.FieldsStatistics.GrainStorage.size();
    vector<double> pointsPerGrain;
    avgCRSS.Allocate({numberOfGrains, 12});
    for (int i = 0; i < numberOfGrains; i++)
    {
        pointsPerGrain.push_back(0.0);
        for(int slipSys = 0; slipSys < 12; slipSys++)
        {
            avgCRSS({i, slipSys}) = 0.0;
        }
    }
    STORAGE_LOOP_BEGIN(i,j,k, Phase.Fields,0)
    {
        if (Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin(); alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                int phaseIndex = alpha->index;
                if (alpha->value > 0.5)
                {
                    pointsPerGrain[phaseIndex] += alpha->value;
                    for(int slipSys = 0; slipSys < 12; slipSys++)
                    {
                        double tmpCRSS =  PFCP.CRSS(i,j,k).get(phaseIndex, slipSys);
                        avgCRSS({phaseIndex, slipSys}) += tmpCRSS * alpha->value;
                    }
                }
            }
        }
        else
        {
            int phaseIndex = Phase.Fields(i,j,k).front().index;
            pointsPerGrain[phaseIndex] += 1.0;
            for(int slipSys = 0; slipSys < 12; slipSys++)
            {
                double tmpCRSS =  PFCP.CRSS(i,j,k).get(phaseIndex, slipSys);
                avgCRSS({phaseIndex, slipSys}) += tmpCRSS;
            }
        }
    }
    STORAGE_LOOP_END

    for (int it = 0; it != numberOfGrains; it++)
    {
        for(int slipSys = 0; slipSys < 12; slipSys++)
        {
            int thPhaseIndex = Phase.FieldsStatistics[it].Phase;
            double phaseCRSS = PFCP.PhaseCriticalResolvedShearStress({thPhaseIndex, slipSys});
            if(pointsPerGrain[it] != 0)
            {
                avgCRSS({it, slipSys}) /= pointsPerGrain[it];
                if (avgCRSS({it, slipSys}) < phaseCRSS)
                {
                    avgCRSS({it, slipSys}) = phaseCRSS;
                }
            }
            else
            {
                avgCRSS({it, slipSys}) = phaseCRSS;
            }
        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields, 0,)
    {
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin(); alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                int phaseIndex = alpha->index;
                int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
                if (alpha->value <= DBL_EPSILON)
                {
                    for (int ss = 0; ss < 12; ss++)
                    {
                        PFCP.CRSS(i,j,k).set(phaseIndex, ss, avgCRSS({phaseIndex, ss}));
                    }
                }
                else
                {
                    for (int ss = 0; ss < 12; ss++)
                    {
                        if (PFCP.CRSS(i,j,k).get(phaseIndex, ss) < PFCP.PhaseCriticalResolvedShearStress({thPhaseIndex, ss}))
                        {
                            PFCP.CRSS(i,j,k).set(phaseIndex, ss, avgCRSS({phaseIndex, ss}));
                        }
                    }
                }
            }
        }
        else
        {
            int phaseIndex = Phase.Fields(i,j,k).front().index;
            int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
            for (int ss = 0; ss < 12; ss++)
            {
                double phaseCRSS = PFCP.PhaseCriticalResolvedShearStress({thPhaseIndex, ss});
                if (PFCP.CRSS(i,j,k).get(phaseIndex, ss) < phaseCRSS)
                {
                    PFCP.CRSS(i,j,k).set(phaseIndex, ss, phaseCRSS);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
PFCP.SetBoundaryConditions(BC);
}
void PlasticFlowMethods::UpdateCRSS(std::vector<OPObject*> Objects)
{
    string fctName = "UpdateCRSS";
    string className = "PlasticFlowMethods";

    PhaseField& Phase           = *(static_cast<PhaseField*>(OPObject::findOPObject(Objects, "PhaseField", className, fctName, true)));
    PlasticFlowCP& PFCP         = *(static_cast<PlasticFlowCP*>(OPObject::findOPObject(Objects, "PlasticFlowCP", className, fctName, true)));
    BoundaryConditions& BC      = *(static_cast<BoundaryConditions*>(OPObject::findOPObject(Objects, "BoundaryConditions", className, fctName, true)));

    STORAGE_LOOP_BEGIN(i,j,k, Phase.Fields,0)
    {
        for(auto alpha = Phase.Fields(i, j, k).cbegin(); alpha < Phase.Fields(i, j, k).cend(); ++alpha)
        {
            if (alpha->value <= DBL_EPSILON )
                for (int ss = 0; ss < 12; ss++)
                {
                    PFCP.CRSS(i,j,k).set(alpha->index, ss, 0.0);
                }
            if (Phase.Interface(i,j,k))
            {
                int phaseIndex = alpha->index;
                int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
                if (alpha->value <= DBL_EPSILON)
                {
                    for (int ss = 0; ss < 12; ss++)
                    {
                        double phaseCRSS = PFCP.PhaseCriticalResolvedShearStress({thPhaseIndex, ss});
                        PFCP.CRSS(i,j,k).set(phaseIndex, ss, phaseCRSS);
                    }
                }
                else
                {
                    for (int ss = 0; ss < 12; ss++)
                    {
                        double phaseCRSS = PFCP.PhaseCriticalResolvedShearStress({thPhaseIndex, ss});
                        if (PFCP.CRSS(i,j,k).get(phaseIndex, ss) < phaseCRSS)
                            PFCP.CRSS(i,j,k).set(phaseIndex, ss, phaseCRSS);
                    }
                }
            }
            else
            {
                int phaseIndex = Phase.Fields(i,j,k).front().index;
                int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
                for (int ss = 0; ss < 12; ss++)
                {
                    double phaseCRSS = PFCP.PhaseCriticalResolvedShearStress({thPhaseIndex, ss});
                    if (PFCP.CRSS(i,j,k).get(phaseIndex, ss) < phaseCRSS)
                        PFCP.CRSS(i,j,k).set(phaseIndex, ss, phaseCRSS);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    PFCP.SetBoundaryConditions(BC);
}

void PlasticFlowMethods::ApplyPhaseTransformationFluxes(std::vector<OPObject*> Objects)
{
    string fctName = "ApplyPhaseTransformationFluxes";
    string className = "PlasticFlowMethods";

    PhaseField& Phase = *(static_cast<PhaseField*>(OPObject::findOPObject(Objects, "PhaseField", className, fctName, true)));
    PlasticFlowCP& PFCP           = *(static_cast<PlasticFlowCP*>(OPObject::findOPObject(Objects, "PlasticFlowCP", className, fctName, true)));
    double dt = PFCP.dt;


    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, Phase.Fields, 0,)
    {
        if(Phase.Fields(i,j,k).flag)
        {
            int numberPF = Phase.FieldsStatistics.GrainStorage.size();
            vector<double> NewFractions(numberPF, 0.0);
            vector<double> OldFractions(numberPF, 0.0);
            for (auto it = Phase.Fields(i,j,k).cbegin();
                    it != Phase.Fields(i,j,k).cend(); ++it)
            {
                NewFractions[it->index] = it->value;
                OldFractions[it->index] = it->value;
            }
            for(auto it = Phase.FieldsDot(i,j,k).begin(); it != Phase.FieldsDot(i,j,k).end(); ++it)
            {
                NewFractions[it->indexA] += it->value*dt;
                NewFractions[it->indexB] -= it->value*dt;
            }
            for(auto it = Phase.FieldsDot(i,j,k).begin(); it != Phase.FieldsDot(i,j,k).end(); ++it)
                for(int ss = 0; ss < 12; ss++)
                {
                        int thPhaseIndexA = Phase.FieldsStatistics[it->indexA].Phase;
                        int thPhaseIndexB = Phase.FieldsStatistics[it->indexB].Phase;
                    double PhaseCRSSalpha = PFCP.PhaseCriticalResolvedShearStress({thPhaseIndexA, ss});
                    double PhaseCRSSbeta  = PFCP.PhaseCriticalResolvedShearStress({thPhaseIndexB, ss});
                     double currCRSSalpha = PFCP.CRSS(i,j,k).get(it->indexA,ss);
                     double currCRSSbeta  = PFCP.CRSS(i,j,k).get(it->indexB,ss);

                   if ( NewFractions[it->indexA] > DBL_EPSILON )
                   {
                       double temp = ( OldFractions[it->indexA] / NewFractions[it->indexA] ) * ( currCRSSalpha - PhaseCRSSalpha);
                       PFCP.CRSS(i,j,k).set(it->indexA, ss, temp + PhaseCRSSalpha);
                   }
                   if ( NewFractions[it->indexB] > DBL_EPSILON )
                   {
                        double temp =  ( OldFractions[it->indexB] / NewFractions[it->indexB] ) * (currCRSSbeta - PhaseCRSSbeta);
                        PFCP.CRSS(i,j,k).set(it->indexB, ss, temp + PhaseCRSSbeta);
                   }
                }
        }
    }OMP_PARALLEL_STORAGE_LOOP_END
}

void PlasticFlowMethods::ApplyPlasticStressFreeStrainContribution(
                                std::vector<OPObject*> Objects, double fraction)
{
    string fctName = "ApplyPlasticStressFreeStrainContribution";
    string className = "PlasticFlowMethods";

    PlasticProperties& PFP  = *(static_cast<PlasticProperties*>(OPObject::findOPObject(Objects, "PlasticProperties", className, fctName, true)));
    ElasticProperties& EP       = *(static_cast<ElasticProperties*>(OPObject::findOPObject(Objects, "ElasticProperties", className, fctName, true)));

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.EffectiveEigenStrains, 0,)
    {
        EP.EffectiveEigenStrains(i,j,k) += PFP.PlasticStrain(i,j,k)*fraction;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

vStrain PlasticFlowMethods::GetAvgIntfPlStrain(PhaseField& Phase, PlasticProperties& PFP, int targetPhaseIndex1, int targetPhaseIndex2)
{
    vStrain avgPlStrain;
    avgPlStrain.set_to_zero();
    double nPoints = 0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, PFP.PlasticStrain, PFP.PlasticStrain.Bcells(), reduction(+:nPoints))
    {
        if (Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin(); alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                if (alpha->index == targetPhaseIndex1)
                {
                    for(auto beta = Phase.Fields(i, j, k).cbegin(); beta < Phase.Fields(i, j, k).cend(); ++beta)
                    {
                        if (targetPhaseIndex2 == beta->index)
                        {
#ifdef _OPENMP
#pragma omp critical
#endif
                            {
                                avgPlStrain += PFP.PlasticStrain(i,j,k);
                            }
                            nPoints++;
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    if (nPoints > 0)
    {
        avgPlStrain /= nPoints;
    }
    return avgPlStrain;
}

vStrain PlasticFlowMethods::GetAvgPhasePlStrain(PhaseField& Phase, PlasticProperties& PFP, int targetPhaseIndex)
{
    vStrain avgPlStrain;
    avgPlStrain.set_to_zero();
    double nPoints = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, PFP.PlasticStrain, PFP.PlasticStrain.Bcells(), reduction(+:nPoints))
    {
        if (Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin(); alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                if (alpha->index == targetPhaseIndex)
                {
#ifdef _OPENMP
#pragma omp critical
#endif
                    {
                        avgPlStrain += PFP.PlasticStrain(i,j,k) * alpha->value;
                    }
                    nPoints += alpha->value;
                }
            }
        }
        else
        {
            if (Phase.Fields(i,j,k).front().index == targetPhaseIndex)
            {
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    avgPlStrain += PFP.PlasticStrain(i,j,k);
                }
                nPoints++;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    if (nPoints > 0)
    {
        avgPlStrain /= nPoints;
    }
    return avgPlStrain;
}

double PlasticFlowMethods::GetAvgIntfCRSS(PhaseField& Phase, PlasticFlowCP& PFCP, int targetPhaseIndex1, int targetPhaseIndex2)
{
    double avgCRSS = 0;;
    double nPoints = 0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, PFCP.CRSS, PFCP.CRSS.Bcells(), reduction(+:avgCRSS) reduction(+:nPoints))
    {
        if (Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin(); alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                if (alpha->index == targetPhaseIndex1)
                {
                    for(auto beta = alpha; beta < Phase.Fields(i, j, k).cend(); ++beta)
                    {
                        if (targetPhaseIndex2 == beta->index)
                        {
                            for (int ss = 0; ss < 12; ss++)
                            {
                                avgCRSS += (PFCP.CRSS(i,j,k).get(alpha->index, ss) + PFCP.CRSS(i,j,k).get(beta->index, ss)) / 12.0;
                            }
                                nPoints++;
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    if (nPoints > 0)
    {
        avgCRSS /= nPoints;
    }
    return avgCRSS;
}

double PlasticFlowMethods::GetAvgIntfCRSSOtherPhase(PhaseField& Phase, PlasticFlowCP& PFCP, int targetPhaseIndex1, int targetPhaseIndex2, int thTargetPhaseIndex)
{
    double avgCRSS = 0;
    double nPoints = 0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, PFCP.CRSS, PFCP.CRSS.Bcells(), reduction(+:avgCRSS) reduction(+:nPoints))
    {
        if (Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin(); alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                if (alpha->index == targetPhaseIndex1)
                {
                    for(auto beta = Phase.Fields(i, j, k).cbegin(); beta < Phase.Fields(i, j, k).cend(); ++beta)
                    {
                        if (targetPhaseIndex2 == beta->index)
                        {
                            for(auto gamma = Phase.Fields(i, j, k).cbegin(); gamma < Phase.Fields(i, j, k).cend(); ++gamma)
                            {
                                int phaseIndex = gamma->index;
                                int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
                                if (thPhaseIndex == thTargetPhaseIndex)
                                {
                                    for (int ss = 0; ss < 12; ss++)
                                    {
                                        avgCRSS += PFCP.CRSS(i,j,k).get(phaseIndex, ss) / 12.0;
                                    }
                                }
                            }
                            nPoints++;
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    if (nPoints > 0)
    {
        avgCRSS /= nPoints;
    }
    return avgCRSS;
}

void PlasticFlowMethods::ApplyPhaseTransformationFluxes(PhaseField& Phase, PlasticFlowCP& PFCP, BoundaryConditions& BC)
{
    string fctName = "ApplyPhaseTransformationFluxes";
    string className = "PlasticFlowMethods";

    double dt = PFCP.dt;
    vector<bool> inherit = {true, true};

    STORAGE_LOOP_BEGIN(i,j,k, Phase.Fields,0)
    if(Phase.Fields(i,j,k).flag)
    {
        int numberPF = Phase.FieldsStatistics.GrainStorage.size();
        vector<double> NewFractions(numberPF, 0.0);
        vector<double> OldFractions(numberPF, 0.0);
        NodeVn<12> NewCRSS;
        for(int n = 0; n < numberPF; n++)
        {
            NewCRSS.set_to_zero(n);
        }
        for (auto it = Phase.Fields(i,j,k).cbegin();
                  it != Phase.Fields(i,j,k).cend(); ++it)
        {
            NewFractions[it->index] = it->value;
            OldFractions[it->index] = it->value;
        }
        for(auto it = Phase.FieldsDot(i,j,k).begin(); it != Phase.FieldsDot(i,j,k).end(); ++it)
        {
            NewFractions[it->indexA] += it->value*dt;
            NewFractions[it->indexB] -= it->value*dt;
        }

        for(auto it = Phase.FieldsDot(i,j,k).begin(); it != Phase.FieldsDot(i,j,k).end(); ++it)
        for(int ss = 0; ss < 12; ss++)
        {
            double CRSSAlpha = PFCP.CRSS(i,j,k).get(it->indexA, ss);
            double CRSSBeta = PFCP.CRSS(i,j,k).get(it->indexB, ss);

            if(    NewFractions[it->indexA] > DBL_EPSILON
                    and NewFractions[it->indexB] > DBL_EPSILON
                    and it->value > 0)
            {
                double temp =
                        it->value*dt*(CRSSBeta - CRSSAlpha)/
                        NewFractions[it->indexA];
                NewCRSS.set(it->indexA, ss, temp);
            }
        }
        NewCRSS = NewCRSS + PFCP.CRSS(i,j,k);
        PFCP.CRSS(i,j,k) = NewCRSS;
    }
    STORAGE_LOOP_END
    PFCP.SetBoundaryConditions(BC);
}

void PlasticFlowMethods::DamageCRSS(PhaseField& Phase, DamageProperties& DP, PlasticFlowCP& PFCP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DP.EffectiveDamage,0,)
    {
        if (Phase.Interface(i,j,k))
        {
            for(auto it = Phase.Fields(i, j, k).cbegin(); it != Phase.Fields(i, j, k).cend(); it++)
            {
                int phaseIndex = it->index;
                dVector<12> temp = PFCP.CRSS(i,j,k).get(phaseIndex);
                temp = temp * (1. - DP.EffectiveDamage(i,j,k));
                PFCP.CRSS(i,j,k).set(phaseIndex, temp);
            }
        }
        else
        {
            int phaseIndex = Phase.Fields(i,j,k).front().index;
            dVector<12> temp = PFCP.CRSS(i,j,k).get(phaseIndex);
            temp = temp * (1. - DP.EffectiveDamage(i,j,k));
            PFCP.CRSS(i,j,k).set(phaseIndex, temp);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
} // namespace openphase
