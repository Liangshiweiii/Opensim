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

#ifndef ELASTICPROPERTIES_H
#define ELASTICPROPERTIES_H

#include "Tools/Includes.h"
#include "SymmetryVariants.h"
namespace opensim
{

class Settings;
class Composition;
class PhaseField;
class Advection;
class BoundaryConditions;
class Temperature;
class Orientations;
class PlasticProperties;
class Velocities;
class DamageProperties;
class DrivingForce;
class Diffusion;
class InterfaceEnergy;

const int FreeBoundaries = 0;
const int AppliedStrain  = 1;
const int AppliedStress  = 2;

const int Khachaturyan = 0;
const int Steinbach    = 1;
const int Voigt        = 2;
const int Reuss        = 3;

class ElasticProperties : public OPObject                                       /// Module which stores and handles elastic properties
{
    //friend class Advection;
 public:
    //~ElasticProperties(void);                                                 ///< Destructor, just to indicate that module did not crash at exit
    ElasticProperties(){};                                                      ///< Constructor (does nothing else)
    ElasticProperties(Settings& locSettings,
            const std::string InputFileName = "default");                       ///< Constructs and initializes class

    void Initialize(Settings& locSettings);                                     ///< Initializes the module, allocate the storage, assign internal variables
    void InitializeLD();
    using OPObject::ReadInput;                                                  ///< Reads elastic properties for each phase
    void ReadInput(std::string InputFileName);                                  ///< Reads elastic properties for each phase

    void SetBoundaryConditions(const BoundaryConditions& BC);                   ///< Sets boundary conditions for strains, stresses and displacements
    void SetGrainsProperties(PhaseField& Phase);                                ///< Sets elastic properties for each grain according to its orientation
    void CalculateThermalStrains(PhaseField& Phase, Temperature& Tx);           ///< Calculates thermal expansion contribution to the eigenstrain

    void SetEffectiveEigenStrains(PhaseField& Phase);                           ///< Sets effective eigenstrains
    void SetEffectiveEigenStrains(PhaseField& Phase, Composition& Cx);          ///< Sets effective eigenstrains considering chemo-mechanical coupling

    void SetEffectiveElasticConstants(PhaseField& Phase);                       ///< Sets effective eigenstrains
    void SetEffectiveElasticConstants(PhaseField& Phase, Composition& Cx);      ///< Sets effective eigenstrains considering chemo-mechanical coupling

    void CalculateDrivingForce(PhaseField& Phase, DrivingForce& dGab);          ///< Calculates elastic driving force
    void CalculateDrivingForceLimited(PhaseField& Phase, DrivingForce& dGab);   ///< Calculates elastic driving force limited by Neuber approximation
    void CalculateChemicalPotentialContribution(PhaseField& Phase,
                                             Diffusion& DF);///< Calculates chemical potential contribution
    void CalculateInterfaceEnergyContribution(PhaseField& Phase, InterfaceEnergy& IE);///< Calculates interface energy contribution
    double Energy(void);                                                        ///< Average elastic energy density in Joule/<simulation cell>
    double PointEnergy(int i, int j, int k);                                    ///< Returns elastic energy in a given point
    double AverageEnergyDensity(void);                                          ///< Average elastic energy density in Joule/metre^3
    void WriteEnergyVTK(int tStep);                                             ///< Writes local elastic energy in VTK format

    void AdvectEigenStrain(Velocities& Vel, BoundaryConditions& BC, double dt,
                                             int scheme = Upwind);              ///< Advects the Eigenstrain
    void AdvectStress(Velocities& Vel, BoundaryConditions& BC, double dt,
                                             int scheme = Upwind);              ///< solves d_t sigma = v*d_x(sigma), incompressible advection of stress
    void AdvectStressCompressible(Velocities& Vel, BoundaryConditions& BC,
                                    double dt, int scheme = Upwind);            ///< solves d_t sigma = d_x(v*sigma), compressible advection of stress

    void ReadStrains(int tStep);                                                ///< Reads strains in binary format
    void ReadStresses(int tStep);                                               ///< Reads stresses in binary format
    void WriteStrains(int tStep);                                               ///< Writes strains in binary format
    void WriteStrainsVTK(const int tStep) const;                                ///< Writes strains in VTK format
    void WriteStrainsVTKData(std::stringstream& buffer) const;                  ///< Writes strains in VTK format into a buffer so that multiple VTK outputs can be put in one file

    void WriteEffectiveEigenStrainsTensorVTK(int tStep);                        ///< Writes EffectiveEigenstrains using tensor index notation, there is a factor 0.5 difference to the Voigt notation
    void WriteEffectiveEigenStrainsVoigtVTK(int tStep);                         ///< Writes EffectiveEigenstrains using Voigt index notation, there is a factor 2 difference to the tensor notation
    void WriteStresses(int tStep);                                              ///< Writes stresses in binary format
    void WriteStressesVTK(const int tStep) const;                               ///< Writes stresses in VTK format
    void WriteStressesVTKData(std::stringstream& buffer) const;                 ///< Writes stresses in VTK format into a buffer so that multiple VTK outputs can be put in one file

    void WriteStressIncrementsVTK(const int tStep) const;                       ///< Writes accumulated stresses in VTK format
    void WriteVelocityGradientVTK(int tStep);                                   ///< Writes velocity gradient to VTK format
    void WriteDeformationGradientVTK(int tStep);                                ///< Writes deformation gradient to VTK format
    void WriteEffectiveElasticConstantsVTK(int tStep);

    void CalculateCompliences(void);                                            ///< Calculate compliences from ElasticConstants
    void PrintPointStatistics(const int x, const int y, const int z) const;     ///< Prints various properties at a given point (x, y, z) to screen
    void Remesh(int newNx, int newNy, int newNz, BoundaryConditions& BC);       ///< Remeshes the elasticity data
    void SmoothEigenstrain();
    void ApplyDamage(DamageProperties& DP);                                     ///< Applies local damage to EffectiveElasticConstants
    void ApplyRotation(Orientations& OR);                                       ///< Applies local rotation to EffectiveElasticConstants

    //moved

    void SetVelocityGradient(double* (&rlDefGrad)[9], double dt);
    void SetRotations(Orientations& OR, double dt);
    void SetStressesRX(PhaseField& Phase, PlasticProperties& PFP,
            BoundaryConditions& BC, std::vector<int> targetPhases);

    // end moved

    ElasticProperties& operator= (const ElasticProperties& rhs);
    vStress WriteStressStrainData(std::string filename, std::string LDflag = "SD");

    vStrain AverageStrain;                                                      ///< Average strain calculated by the spectral solver.
    vStrain StrainToRemesh;                                                     ///< Accumulates the strain that has to be applied by remeshing
    vStrain RemeshedStrain;                                                     ///< Accumulates the strain that has been applied by remeshing
    dVector6 AppliedStrainBool;                                                 ///< BC flag for coupled boundary conditions (0 if applied strain given in particular component)
    vStress AppliedStress;                                                      ///< Applied stress tensor
    vStrain AppliedStrain;                                                      ///< Applied strain tensor
	vStrain AppliedStrainRate;                                                  ///< Applied strain rate tensor
    vStress EffectiveAppliedStress;
    vStrain EffectiveAppliedStrain;                                             ///< In the small strain limit: applied strain tensor; in the large deformations: incremenent of the applied strain
    vStrain AppliedStrainOLD;                                                   ///< Applied strain tensor from the previous time step
    vStrain AppliedStressOLD;                                                   ///< Applied strain tensor from the previous time step
    //Elasticity Tensors:
    dMatrix6x6   AverageElasticConstants;                                       ///< Average elastic constants of the whole system
    dMatrix6x6   MAXElasticConstants;                                           ///< Maximum elastic constants of the system
    dMatrix6x6   AverageCompliences;                                            ///< Inverse of the AverageElasticConstants
    dMatrix6x6   MINCompliences;                                                ///< Inverse of the MAXElasticConstants
    dVector6     AvgStrainMask;                                                 ///< Diagonal Components are 1 if Free Boundaries is selected for the corresponding direction
    dVector6     LoadStressMask;                                                ///< Diagonal Components are 1 if Applied Stress is selected for the corresponding direction
    dVector6     AppStrainMask;

    Storage<vStrain>       EigenStrains;                                        ///< Reference eigenstrains for each phase field
    Storage<dMatrix6x6>    ElasticConstants;                                    ///< Reference elastic constants for each phase field
    Storage<dMatrix6x6>    Compliences;                                         ///< Inverse of the ElasticConstants

    Storage<vStrain>       PhaseEigenStrains;                                   ///< Reference eigenstrains of each phase (relative to Cref), grain orientation not considered
    Storage<dMatrix6x6>    PhaseElasticConstants;                               ///< Reference elastic constants of each phase (relative to Cref), grain orientation not considered
    Storage<dMatrix6x6>    PhaseCompliences;                                    ///< Inverse of the PhaseElasticConstants, grain orientation not considered
    bool ThermoMechanicalCoupling;                                              ///< Thermo-mechanical coupling flag
	bool NeuberCorrection;														///< Neuber correction flag
    Storage<vStrain>        PhaseAlpha;                                         ///< Thermal expansion coefficients of each phase
    Storage<vStrain>        Alpha;                                              ///< Thermal expansion coefficients of each phase field
    Storage<double>         Tref;                                               ///< Thermal expansion reference temperature of each phase


    /* template specialization integer number stands for
     * the number of extra dimensions
     * order:
     * Rank = 2: phase, component
     * Rank = 1: phase or component
     */

    bool ChemoMechanicalCoupling;                                               ///< Chemo-mechanical coupling flag
    std::vector < std::string > Names;                                          ///< Names of corresponding chemical components (e.g Au, Cu, Na, Cl etc.)

    Tensor<dMatrix6x6, 2>   Kappa;                                              ///< Coupling constants for the elasticity and diffusion for each phase field
    Tensor<dMatrix6x6, 2>   PhaseKappa;                                         ///< Coupling constants for the elasticity and diffusion for each phase
    Tensor<vStrain, 2>      Lambda;                                             ///< Coupling constants for the eigenstrain and diffusion for each phase field
    Tensor<vStrain, 2>      PhaseLambda;                                        ///< Coupling constants for the eigenstrain and diffusion for each phase
    Tensor<double, 2>       Cref;                                               ///< Reference compositions for each phase

    Storage3D<dMatrix3x3, 0>    VelocityGradient;                               ///< Storage for deformation gradient rate matrix
    Storage3D<vStress, 0>       Stresses;                                       ///< Storage for stresses
    Storage3D<vStress, 0>       StressIncrements;                               ///< Storage for stress increments
    Storage3D<vStrain, 0>       Strains;                                        ///< Storage for strains
    Storage3D<vStrain, 0>       StrainIncrements;                               ///< Storage for strain increments

    Storage3D<vStrain, 0>       EffectiveEigenStrains;                          ///< Storage for effective eigenstrains
    Storage3D<dMatrix6x6, 0>    EffectiveElasticConstants;                      ///< Storage for effective elastic constants

    Storage3D<vStrain, 0>       StrainsDot;                                     ///< Storage for the increment of the strains due to advection
    Storage3D<vStress, 0>       StressesDot;                                    ///< storage for the increment of the stresses due to advection
    SymmetryVariants            Variants;                                       ///< Symmetry variants storage

    int Nx;
    int Ny;
    int Nz;
    double dx;
	double dt;
    int Nphases;
    int Ncomp;

    int ElasticityModel;                                                        ///< Elasticity model selector

    double DDInitial;
    double DDSaturation;
    
    bool KeepAspectRatio;                                                       ///< Restricts Free Boundary conditions to only volume relaxation
    bool LargeDeformations;                                                     ///< Large deformation flag.

 protected:
 private:

};
}// namespace openphase
#endif
