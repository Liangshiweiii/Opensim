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

#ifndef HEATDIFFUSION_H
#define HEATDIFFUSION_H

#include "Tools/Includes.h"

namespace opensim
{
class Settings;
class GrainInfo;
class PhaseField;
class Temperature;
class BoundaryConditions;

class HeatDiffusion : public OPObject                                           ///< Heat equation solver class
{
 public:
    HeatDiffusion(){};                                                          ///< Constructor
    HeatDiffusion(const Settings& locSettings);                                 ///< Constructor which directly initializes class from input file
    HeatDiffusion(const Settings& locSettings, const std::string FileName);     ///< Constructor which directly initializes class from input file

    using OPObject::ReadInput;
    void Initialize(const Settings& locSettings);                               ///< Allocates memory, initializes the settings);
    void Initialize(const Settings& locSettings, const std::string FileName);   ///< Allocates memory, initializes the settings, reads input
    void ReadInput(const std::string FileName);                                 ///< Reads input parameters
    void Remesh(const int newNx, const int newNy, const int newNz);             ///< Change system size while keeping the data
    void Solve(const PhaseField& Phase, Temperature& Temperature,
            const BoundaryConditions& BC, const double dt);                     ///< Solver

    double MaxThermalDiffusivity;                                               ///< Maximal thermal diffusivity in the system

    Storage<double> PhaseThermalDiffusivity;                                    ///< Thermal diffusivity
    Storage<double> PhaseHeatCapacity;                                          ///< Heat capacity
    //Storage<double> PhaseDensity;                                             ///< NOT USED! Density

    Storage3D<double,0> EffectiveThermalDiffusivity;                            ///< Thermal diffusivity
    Storage3D<double,0> EffectiveHeatCapacity;                                  ///< Heat capacity
    //Storage3D<double,0> EffectiveDensity;                                     ///< NOT USED! Density

 protected:
 private:

    int Nx;                                                                     ///< Size of the inner calculation domain along X
    int Ny;                                                                     ///< Size of the inner calculation domain along Y
    int Nz;                                                                     ///< Size of the inner calculation domain along Y
    int Nphases;                                                                ///< Number of phases
    double dx;                                                                  ///< Grid spacing
    double dx2Inv;

    void CalculateLaplacian(Temperature& Temperature) const;                    ///< Calculation of temperature gradients using finite differences
    void ClearEffectiveProperties();                                            ///< Clears effective properties storages
    void SetEffectiveProperties(const PhaseField& Phase);                       ///< Calculation of temperature gradients using finite differences
    void UpdateTemperature(Temperature& Temperature, const double dt) const;    ///< Explicit time integration of temperature
};

} // namespace opensim
#endif
