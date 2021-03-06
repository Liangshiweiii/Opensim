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

#ifndef INTERFACEDIFFUSION_H
#define INTERFACEDIFFUSION_H

#include "Tools/Node.h"
#include "Tools/NodeV.h"

namespace opensim
{
class InterfaceEnergy;
class PhaseField;
class Settings;

class InterfaceDiffusion : public OPObject                                      ///< Calculates the interface diffusion inside the interface of multiple phases
{
 public:
    InterfaceDiffusion(){};                                                     ///< Empty constructor
    InterfaceDiffusion(const Settings& locSetting);                             ///< Constructor that initialises the class form input file

    using OPObject::Initialize;                                                 ///< Rads the initial parameters form default path/filename
    void Initialize(const Settings& locSettings);                               ///< Initializes, just to indicate that module has been created.

    using OPObject::ReadInput;                                                  ///< Rads the initial parameters form default path/filename
    void ReadInput(const std::string InputFileName);                            ///< Reads the initial parameters for the interface diffusion

    void CalculatePhaseFieldIncrements(PhaseField& Phase,
            const InterfaceEnergy& Sigma);                                      ///< Calculates the movement of the interface due to the surface diffusion

    double PotentialDerivative(const double alpha, const double beta) const;
    void CalculateDiffusionPotential(PhaseField& Phase, const InterfaceEnergy& Sigma);///< Calculates diffusion potential to the interface curvature
    void CalculateDiffusionPotentialGradients();                                ///< Calculates the gradient/diffusion-flux
    void CalculateDiffusionPotentialLaplacian(PhaseField& Phase);               ///< Laplacian of the diffusion potential

    Matrix<double>        Coefficients;                                         ///< Interface diffusion coefficient of each interface
    Storage3D< Node,  0 > DiffusionPotential;                                   ///< Interface diffusion-potential of the tangential diffusion-flux
    double                dx;                                                   ///< Grid spacing x-direction of the system
    double                dy;                                                   ///< Grid spacing x-direction of the system
    double                dz;                                                   ///< Grid spacing x-direction of the system
    double                DoubleObstacleSmoothnessRange;                        ///< Defines the smoothing of the double obstacle potential
    unsigned int          Nphases;                                              ///< Number of thermodynamic phases
    unsigned int          Nx;                                                   ///< X dimension of the system
    unsigned int          Ny;                                                   ///< Y dimension of the system
    unsigned int          Nz;                                                   ///< Z dimension of the system

 protected:
 private:
};

class InterfaceDiffusionAnisotropic : public InterfaceDiffusion                 ///< Calculates the interface diffusion inside the interface of multiple phases
{
 public:
    InterfaceDiffusionAnisotropic(){};                                          ///< Empty constructor
    InterfaceDiffusionAnisotropic(const Settings& locSetting);                  ///< Constructor that initialises the class form input file
    void Initialize(const Settings& locSettings);                               ///< Initializes, just to indicate that module has been created.
    void ReadInput(const std::string InputFileName);                            ///< Reads the initial parameters for the interface diffusion

    void CalculatePhaseFieldIncrements(PhaseField& Phase,
            const InterfaceEnergy& Sigma);                                      ///< Calculates the movement of the interface due to the surface diffusion
 protected:
    double PotentialDerivative(const double alpha, const double beta) const;
    void CalculateDiffusionFlux(PhaseField& Phase);                             ///< Calculates the interface diffusion flux
    void CalculateDiffusionFluxDivergence(PhaseField& Phase);                   ///< Calculates the divergence of the diffusion potential
    void CalculateDiffusionPotentialGradients();                                ///< Calculates the gradient/diffusion-flux

    Storage3D< NodeV, 0 > DiffusionFlux;                                        ///< Interface diffusion-potential of the tangential diffusion-flux
 private:
};
}
#endif

