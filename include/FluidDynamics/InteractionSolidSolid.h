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

#ifndef INTERACTIONSOLIDSOLID_H
#define INTERACTIONSOLIDSOLID_H

#include "PhaseField.h"

namespace openphase
{

class D3Q27;
class FlowSolverLBM;
class Settings;
class Velocities;

class InteractionSolidSolid : public OPObject                                   ///<  Handles the interaction between solid particles in the fluid environment
{
 public:

    using OPObject::ReadInput;
    void Initialize(Settings& locSettings);                                     ///<  Allocates memory and initiates global settings, structure described in OPSettings.h file
    void ReadInput(std::string InputFileName);                                  ///<  Reads input parameters from a file

    static void CalculateSolidSolid(PhaseField& Phase,
            const int i, const int j, const int k,
            const double dt, const int order = 4, const int cutoff = 3,
            const double strength = 1., const double elastic = 0.);             ///< Calculates the solid-solid interaction at the point (i,j,k)

    static void CalculateSolidSolid(PhaseField& Phase,
            const FlowSolverLBM& LBM,
            const double dt, const int order = 4, const int cutoff = 3,
            const double strength = 1., const double elastic = 0.);             ///< Calculates the solid-solid interaction the entire simulation domain, but only near LMB.Obstacle!

    static void CalculateSolidSolid(PhaseField& Phase,
            const double dt, const int order = 4, const int cutoff = 3,
            const double strength = 1., const double elastic = 0.);             ///< Calculates the solid-solid interaction the entire simulation domain

    static void AdvectSolid(PhaseField& Phase, const BoundaryConditions& BC,    ///< NOTE: deprecated use AdvectionHR!
            const double dt);

    static void PreserveVolume(PhaseField& Phase,
            const std::vector<double>& RefVolume, const int _index,
            const double dt);                                                   ///< NOTE: deprecated use AdvectionHR!

    static void SetRefVolume(const PhaseField& Phase,
            std::vector<double>& RefVolume);                                    ///< NOTE: deprecated use AdvectionHR!

    static double VolumeError(const PhaseField& Phase,
            const std::vector<double>& RefVolume, const int index);             ///< NOTE: deprecated use AdvectionHR!
 protected:
 private:
};


} //namespace openphase
#endif
