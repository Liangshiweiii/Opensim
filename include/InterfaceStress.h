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

#ifndef INTERFACESTRESS_H
#define INTERFACESTRESS_H

#include "Tools/OPObjects.h"
namespace opensim
{
class ElasticProperties;
class InterfaceEnergy;
class PhaseField;
class vStress;

struct InterfaceStress : public OPObject                                        ///<  Handles isotropic interface stresses for each interface
{
    InterfaceStress(const Settings& locSetting);                                ///< Constructor that initialises the class form input file
    using OPObject::Initialize;
    void Initialize(const Settings& locSettings);                               ///< Initializes the module, allocate the storage, assign internal variables

    const vStress operator()(const long int i, const long int j,
            const long int k, const PhaseField& Phase,
            const InterfaceEnergy& Sigma) const;                                ///< Calculates interface stress at the point (i,j,k)

    void SetEffectiveEigenStrains(ElasticProperties& EP,
            const PhaseField& Phase,
            const InterfaceEnergy& Sigma) const;                                ///< Calculates the resulting volume strain due to the stress of the interface
};

}// namespace openphase
#endif
