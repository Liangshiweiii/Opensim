/* This file is part of OPENPHASE BASE.
 * Copyright (c) 2020 OpenPhase Solutions GmbH, Bochum, Germany
 * For more details visit https://www.openphase-solutions.com
 * 
 * OPENPHASE BASE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *  
 * OPENPHASE BASE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *  
 * You should have received a copy of the GNU General Public License
 * along with OPENPHASE BASE.  If not, see <http://www.gnu.org/licenses/>. */

#ifndef LARGEDEFORMATIONSMETHODS_H
#define LARGEDEFORMATIONSMETHODS_H

#include "fftw3.h"
#include "Tools/Includes.h"
#include "Tools/NodeVn.h"

namespace opensim
{
class ElasticProperties;
class Orientations;
class Velocities;
class Settings;
class DrivingForce;
class SpectralElasticSolver;
class PhaseField;
class Advection;
class Composition;
class InterfaceField;
class BoundaryConditions;
class InterfaceEnergy;
class InterfaceMobility;
class PlasticProperties;
class DamageProperties;
class PlasticFlowCP;

class LargeDeformationsMethods
{
 public:
    static void AccumulateStresses(ElasticProperties& EP, double dt);           /// Adds stress increment coming from the elastic solver to the accumulated stress using the Truesdell rate formula
    static void AccumulateStressesRotation(ElasticProperties& EP,
                       Orientations& OR, Storage3D<dMatrix3x3, 0>& RotationsOLD,
                         /*Storage3D<dMatrix6x6, 0>& StiffnessOLD,*/ double dt);
    static void RotatePlasticStrains(PlasticProperties& PFP,
					   Orientations& OR, Storage3D<dMatrix3x3, 0>& RotationsOLD,
						 /*Storage3D<dMatrix6x6, 0>& StiffnessOLD,*/ double dt);
    static void ApplyHomExpansionToAccStress(ElasticProperties& EP,
                                                             vStrain StrainInc);/// Corrects the accumulated stress according to volume change of the system
    //static void CalculateAverageVelocity(ElasticProperties& EP, Velocities& Vel,
    //                              double dt);                                   /// Calculates the velocity from the displacement increment
    static void CalculateVelocityGradient(ElasticProperties& EP,
                                                               Velocities& Vel);/// Calculates the velocity gradient tensor L
    //static void CalculateDeformationGradient(ElasticProperties& EP);
//    void CalculateAngularVelocity(ElasticProperties& EP);                     /// Calculates the angular velocity from the velocity gradient
    static void SetHomogeneousVelocity(PhaseField& Phase, Velocities& Vel,
                                          vStrain HomogeneousStrain, double dt);
    static void SetHomogeneousVelocityBC(Velocities& Vel);
    static void SaveRotations(Orientations& OR, BoundaryConditions& BC,
                                        Storage3D<dMatrix3x3, 0>& RotationsOLD);
    static void SaveEffectiveEigenStrains(PhaseField& Phase,
                                 ElasticProperties& EP, BoundaryConditions& BC,
                               Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD);/// calculates effective eigenstrain and stores it in EP.EffectiveEigenStrainsOLD
    static void SaveEffectiveEigenStrains(PhaseField& Phase,
                ElasticProperties& EP, Orientations& OR, BoundaryConditions& BC,
                               Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD);
    static void SaveEffectiveEigenStrains(PhaseField& Phase, Composition& Cx,
                                 ElasticProperties& EP, BoundaryConditions& BC,
                               Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD);
    //static void SavePlasticStrain(PlasticProperties& PFP,
    //                                   Storage3D<vStrain, 0>& PlasticStrainOLD);
    static void SaveCRSS(PlasticFlowCP& PFCP, Storage3D<NodeVn<12>, 0>& CRSSOLD);
    static void SaveDamage(DamageProperties& DP,
                                               Storage3D<double, 0>& DamageOLD);
	static void SaveDamageN(DamageProperties& DP,
		Storage3D<double, 0>& DamageOLD, Storage3D<double, 0>& DamageNOLD);
    static void SubtractPlasticStrainOLD(PlasticProperties& PFP,
                                       Storage3D<vStrain, 0>& PlasticStrainOLD);
    static void ADDPlasticStrainOLD(PlasticProperties& PFP,
                                       Storage3D<vStrain, 0>& PlasticStrainOLD);
    static void DeletePlasticStrain(PlasticProperties& PFP);
    //static void RestorePlasticStrain(PlasticProperties& PFP,
    //                                   Storage3D<vStrain, 0>& PlasticStrainOLD);
    static void RestoreCRSS(PlasticFlowCP& PFCP,
    										Storage3D<NodeVn<12>, 0>& CRSSOLD);
    static void RestoreDamage(DamageProperties& DP,
    											Storage3D<double, 0>& DamageOLD);
	static void RestoreDamageN(DamageProperties& DP,
		Storage3D<double, 0>& DamageOLD, Storage3D<double, 0>& DamageNOLD);
	static void RestoreRotations(Orientations& OR,
		BoundaryConditions& BC, Storage3D<dMatrix3x3, 0>& RotationsOLD);
    static void SaveEffectiveKirkendallEigenStrains(PhaseField& Phase,
                Composition& Cx, ElasticProperties& EP, BoundaryConditions& BC,
                               Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD);
    static void SaveEffectiveElasticConstants(ElasticProperties& EP,
                                        Storage3D<dMatrix6x6, 0>& StiffnessOLD);
    static void SavePhaseField(PhaseField& Phase, PhaseField& PhaseOLD);
//    static double GetEigenstrainNorm(PhaseField& Phase, ElasticProperties& EP,
//                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD,
//                                              double MaxAllowedStrainIncrement);
//    static double GetEigenstrainNorm(PhaseField& Phase, ElasticProperties& EP,
//                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD,
//                                CreepProperties& CP, CreepModel& CM,
//                                              double dt, double MaxAllowedStrainIncrement);
    static double GetEigenstrainNorm(std::vector<OPObject*> Objects,
                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD,
                                         Storage3D<dMatrix6x6, 0>& StiffnessOLD,
                                   double dt, int nSubsteps, double MaxAllowedStrainIncrement,
								   bool verbose = true, bool screenWrite = true);
//    static double GetEigenstrainNorm(PhaseField& Phase, ElasticProperties& EP,
//                                                               Orientations& OR,
//                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD,
//                                              double MaxAllowedStrainIncrement);
    static double GetEigenstrainNorm(PhaseField& Phase, Composition& Cx,
                                                          ElasticProperties& EP,
                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD,
                                              double MaxAllowedStrainIncrement);

//    static double SetEigenstrainsIncrement(PhaseField& Phase,
//         ElasticProperties& EP, Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD,
//                                                                  int substeps);
//    static double SetEigenstrainsIncrement(PhaseField& Phase,
//         ElasticProperties& EP, Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD,
//         CreepProperties& CP, CreepModel& CM, double dt, int substeps);
    static double SetEigenstrainsIncrement(std::vector<OPObject*> Objects,
                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD,
                                         Storage3D<dMatrix6x6, 0>& StiffnessOLD,
                                                       double dt, int nSubsteps,
													   bool verbose = true);
//    static double SetEigenstrainsIncrement(PhaseField& Phase,
//                                        ElasticProperties& EP, Orientations& OR,
//                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD,
//                                                                  int substeps);
    static double SetEigenstrainsIncrement(PhaseField& Phase, Composition& Cx,
                                                          ElasticProperties& EP,
                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD,
                                                                  int substeps);
    static double GetAppliedStrainNorm(ElasticProperties& EP,
                                              double MaxAllowedStrainIncrement,
											  bool screenWrite);
    static double GetAppliedStressNorm(ElasticProperties& EP,
                                              double MaxAllowedStrainIncrement);
    static double SetAppliedStrainIncrement(ElasticProperties& EP,
                                                                  int substeps);
    static double SetAppliedStressIncrement(ElasticProperties& EP,
                                                                  int substeps);
    static void SetEffectiveElasticConstants(PhaseField& Phase, ElasticProperties& EP, Orientations& OR);// corrects stiffness tensor according to: Wooseok Ji, Anthony M. Waas, Zdenek P. Bazant - Journal of Applied Mechanics 77(2010) - Errors Caused by Non-Work-Conjugate Stress and Strain Measures and Necessary Corrections in Finite Element Programs
//    double SetKirkendallEigenstrainsIncrement(PhaseField& Phase, Composition& Cx,
//                                      ElasticProperties& EP, int substeps = -1);
    static void RestoreInterface(PhaseField& Phase, BoundaryConditions& BC,
                                  InterfaceEnergy& Sigma, InterfaceMobility& Mu,
                                                InterfaceField& Psi, double dt);// apllies double obstacle potential at the outer interface cells to suppress interface spreading caused by advection and remeshing
    static void ConsiderStiffnessIncrement(ElasticProperties& EP,
                                        Storage3D<dMatrix6x6, 0>& StiffnessOLD);
    static void SetPlasticStrainEquivalent(
                                  Storage3D<double, 0>& PlasticStrainEquivalent,
                                          Storage3D<vStrain, 0>& PlasticStrain);
    static void SetPlasticStrainEquivalentIntf(PhaseField& Phase,
                                  Storage3D<double, 0>& PlasticStrainEquivalent,
                          Storage3D<vStrain, 0>& PlasticStrain, int phaseIndex);
    static void CalculatePhaseDamageEnergy(
                   std::vector<double>& PhaseDamageEnergy, DamageProperties& DP,
                                       ElasticProperties& EP, PhaseField& Phase,
                                               Storage3D<double, 0>& DamageOLD);
    static void WriteEffectiveElasticConstantsDiffVTK(
    										   ElasticProperties& EP, int tStep,
										  Storage3D<dMatrix6x6,0>& StiffnessOld);

 private:
};
} // namespace openphase
#endif //LARGEDEFORMATIONSMETHODS_H
