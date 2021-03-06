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

#ifndef ELASTICITYSOLVER_H
#define ELASTICITYSOLVER_H

#include "Tools/Includes.h"

namespace opensim
{
class ElasticProperties;
class Orientations;
class Settings;
class BoundaryConditions;

class ElasticitySolver : public OPObject                                        ///< Elastic problem solver base class
{
 public:
    virtual ~ElasticitySolver(void){};                                          ///< Destructor
    virtual void Initialize(Settings& locSettings) = 0;                         ///< Named constructor
    virtual void ReInitialize(ElasticProperties& EP) = 0;                       ///< Needed if the system has been remeshed

    virtual int Solve(ElasticProperties& EP, Orientations& OR, BoundaryConditions& BC,
              double StrainsAccuracy, double PressureAllowed, int MAXIterations,
              double dt, bool getU = false) = 0;                                                   ///< Solves elastic problem
    virtual int Solve(ElasticProperties& EP, BoundaryConditions& BC,
                  double StrainsAccuracy, double PressureAllowed, int MAXIterations,
                  double dt, bool getU = false) = 0;
 private:
};
} // namespace openphase
#endif //ELASTICITYSOLVER
