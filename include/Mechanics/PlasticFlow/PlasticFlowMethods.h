

#ifndef PLASTICFLOWMETHODS_H
#define PLASTICFLOWMETHODS_H

#include "Tools/Includes.h"

namespace opensim
{
class Settings;
class PhaseField;
class ElasticProperties;
class Composition;
class Temperature;
class Settings;
class Orientations;
class BoundaryConditions;
class PlasticFlowCP;
class PlasticProperties;
class Composition;
class DamageProperties;

class PlasticFlowMethods                                                        /// Class containing calculation methods for the PlasticFlow model class
{
public:
    static double CalculatePlasticStrainCP(std::vector<OPObject*> Objects, bool verbose);
    static void ApplyHardening(PhaseField& Phase, PlasticFlowCP& PF,
            PlasticProperties& PFP, dVector<12> shearRate,
            int phaseIndex, int i, int j, int k);
    static void ApplyHardening(PhaseField& Phase, PlasticFlowCP& PF,
            PlasticProperties& PFP, Composition& Cx, dVector<12> shearRate,
            int phaseIndex, int i, int j, int k);
    static void UpdateCRSS(const PhaseField& Phase, PlasticFlowCP& PFCP, BoundaryConditions& BC);
    static void UpdateCRSS(std::vector<OPObject*> Objects);
    static void ApplyPhaseTransformationFluxes(std::vector<OPObject*> Objects);
    static void ApplyPlasticStressFreeStrainContribution(
                           std::vector<OPObject*> Objects, double fraction=1.0);
    static vStrain GetAvgIntfPlStrain(PhaseField& Phase, PlasticProperties& PFP, int targetPhaseIndex1, int targetPhaseIndex2);
    static vStrain GetAvgPhasePlStrain(PhaseField& Phase, PlasticProperties& PFP, int targetPhaseIndex);
    static double GetAvgIntfCRSS(PhaseField& Phase, PlasticFlowCP& PFCP, int targetPhaseIndex1, int targetPhaseIndex2);
    static double GetAvgIntfCRSSOtherPhase(PhaseField& Phase, PlasticFlowCP& PFCP, int targetPhaseIndex1, int targetPhaseIndex2, int thTargetPhaseIndex);
    static void ApplyPhaseTransformationFluxes(PhaseField& Phase, PlasticFlowCP& PFCP, BoundaryConditions& BC);

    static void DamageCRSS(PhaseField& Phase, DamageProperties& DP, PlasticFlowCP& PFCP);
};
}// namespace openphase

#endif
