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

#ifndef LARGEDEFORMATIONS_H
#define LARGEDEFORMATIONS_H

#include "fftw3.h"
#include "Tools/Includes.h"
#include "Tools/NodeVn.h"
#include "BoundaryConditions.h"

namespace opensim
{
class ElasticProperties;
class Orientations;
class Velocities;
class Settings;
class DrivingForce;
class PlasticProperties;
class PlasticityModel;
class ElasticSolver;
class PhaseField;
class AdvectionPhaseField;
class AdvectionElasticProperties;
class Composition;
class Temperature;
class BoundaryConditions;

class LargeDeformations : public OPObject
{
 public:
    LargeDeformations(){};
    LargeDeformations(Settings& locSettings, ElasticProperties& EP);
    ~LargeDeformations(void);

    void Initialize(Settings& locSettings, ElasticProperties& EP);
    void ReadInput(std::string FileName);
	using OPObject::ReadInput;
    void SetBoundaryConditions(const BoundaryConditions& BC);
    void Remesh(int Nx, int Ny, int Nz, BoundaryConditions& BC);
    double dx;

    void Solve(std::vector<OPObject*> Objects, bool screenWrite);
	int calculateLDStep(std::vector<OPObject*> Objects, int nSubsteps, bool force, bool screenWrite);

//    std::vector<double> PhaseDamageEnergy;
//    std::vector<double> PhasePlasticEnergy;

 private:
    int counter;
    int maxIterations;
    int minIterations;
	double  dt;

    bool Verbose;
    double MaxAllowedStrainIncrement;
    double StrainAccuracy;
    double MaxPressure;
    int MaxSolverIterations;
    bool IntfStabilisationNeeded;                                               /// Controls whether interface stabilisation is done or not. In case there is no phase evolution calculation advection and remeshing lead to interface spreading. If phase evolution is calculated the problem does not appear and interface stabilisation might lead to stronger grain shrinkage.
    void CalculateAverageVelocity(ElasticProperties& EP, Velocities& Vel, double dt); /// Calculates velocity from the displacement increment
    
    void SavePlasticStrain(PlasticProperties& PFP,
                                       Storage3D<vStrain, 0>& PlasticStrainOLD);
    void RestorePlasticStrain(PlasticProperties& PFP,
                                       Storage3D<vStrain, 0>& PlasticStrainOLD);
    BoundaryConditions VelocityBC;

    Storage3D<vStrain, 0>       EffectiveEigenStrainsOLD;                       /// Storage for effective eigenstrains
    Storage3D<vStrain, 0>       PlasticStrainOLD;                               /// Storage for plastic strain
	Storage3D<vStrain, 0>       PlasticStrainDamage;                            /// Storage for damage calculation
    Storage3D<dMatrix3x3, 0>    RotationsOLD;                                   /// Storage for Rotations
    Storage3D<dMatrix6x6, 0>    StiffnessOLD;                                   /// Storage for ElasticConstants from previous time step
    Storage3D<double, 0>        DamageOLD;
	Storage3D<double, 0>        DamageNOLD;
    Storage3D<NodeVn<12>, 0> 	CRSSOLD;
};
} // namespace openphase
#endif //LARGEDEFORMATIONS_H
