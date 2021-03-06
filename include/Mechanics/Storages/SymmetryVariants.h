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

#ifndef SYMMETRYVARIANTS_H
#define SYMMETRYVARIANTS_H

#include "Tools/Includes.h"
#include "Tools/ElasticityTensors.h"

namespace opensim
{

class Settings;

class SymmetryVariants : OPObject                                               ///< Stores transformation matrices for symmetry variants of thermodynamic phases.
{
 public:
    SymmetryVariants(){};
    SymmetryVariants(const Settings& locSettings, 
                     std::string InputFileName = "default");                    ///< Constructor, uses Initialize() and ReadInput()
    void Initialize(const Settings& locSettings);                               ///< Initializes storages, sets internal variables.
    using OPObject::ReadInput;                                                  ///< Reads the input from the default input file
    void ReadInput(std::string InputFileName);                                  ///< Reads the input from the file InputFileName
    dMatrix3x3& operator()(int PhaseIndex, int VariantIndex)                    ///< Bi-directional access operator
    {
        return TransformationMatrices[PhaseIndex][VariantIndex];
    }
    unsigned int Nphases;                                                       ///< Number of thermodynamic phases
    bool set;                                                                   ///< Indicates if symmetry variants were set
    int Nvariants(const int PhaseIndex) const                                   ///< Returns number of symmetry variants of a given phase
    {
        return TransformationMatrices[PhaseIndex].size();
    }
 protected:
     std::vector< std::vector < dMatrix3x3 > > TransformationMatrices;          ///< Transformation matrices storage
 private:
};

} // namespace openphase

#endif
