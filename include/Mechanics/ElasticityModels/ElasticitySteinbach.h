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

#ifndef ELASTICITYSTEINBACH_H
#define ELASTICITYSTEINBACH_H

namespace opensim
{
class Composition;
class DamageProperties;
class DrivingForce;
class ElasticProperties;
class EquilibriumPartitionDiffusion;
class Orientations;
class PhaseField;
class Plasticity;
class Settings;


class ElasticitySteinbach                                                       ///< Steinbach's approximation for elasticity parameters
{
 public:

    /// Calculates and sets effective eigenstrain in ElasticProperties
    static void SetEffectiveEigenStrains(PhaseField& Phase, ElasticProperties& EP);
    static void SetEffectiveEigenStrains(PhaseField& Phase, ElasticProperties& EP, Orientations& OR);
    static void SetEffectiveEigenStrains(PhaseField& Phase, ElasticProperties& EP, Composition& Cx);
    static void SetEffectiveEigenStrains(PhaseField& Phase, ElasticProperties& EP, Composition& Cx, Orientations& OR);
    /// Calculates and sets effective elastic constants in ElasticProperties
    static void SetEffectiveElasticConstants(PhaseField& Phase, ElasticProperties& EP);
    static void SetEffectiveElasticConstants(PhaseField& Phase, ElasticProperties& EP, Orientations& OR);
    static void SetEffectiveElasticConstants(PhaseField& Phase, ElasticProperties& EP, Orientations& OR, DamageProperties& DP);
    static void SetEffectiveElasticConstants(PhaseField& Phase, ElasticProperties& EP, Composition& Cx);

    static void CalculateDrivingForce(PhaseField& Phase, ElasticProperties& EP, DrivingForce& dGab);
    
    static void CalculateChemicalPotentialContribution(PhaseField& Phase,
                                           ElasticProperties& EP,
                                           EquilibriumPartitionDiffusion& DF);
    //static void CalculateChemicalPotentialContrib(PhaseField& Phase,
    //                          ElasticProperties& EP, ThermodynamicProperties& TP);

 protected:
 private:
};
}// namespace opensim
#endif
