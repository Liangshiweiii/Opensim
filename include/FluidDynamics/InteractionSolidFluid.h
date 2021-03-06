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

#ifndef INTERACTIONSOLIDFLUID_H
#define INTERACTIONSOLIDFLUID_H

#include "Base/Includes.h"
namespace openphase
{

class Settings;
class PhaseField;
class Velocities;
class BoundaryConditions;
class PhaseField;
class InteractionSolidFluid                                                     /// Processes solid liquid interaction
{
 public:
    static void CalculateSolidVelocities(PhaseField& Phase, Velocities& Vel, double dt); /// Calculates velocities of solid particles
    static void CollectGrainsStatistics(PhaseField& Phase);                     ///<  Collects center of mass position and moments of inertia for each phase field.
    static void CollectGrainsStatisticsStepOne(PhaseField& Phase);              ///<  Substep of CollectGrainStatistics
    static void CollectGrainsStatisticsStepTwo(PhaseField& Phase);              ///<  Substep of CollectGrainStatistics
    static void CalculateCenterOfMassWithPeriodicBoundaryConditions(PhaseField& Phase);

 protected:

 private:
};

} //namespace openphase
#endif
