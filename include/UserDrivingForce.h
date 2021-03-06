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

#ifndef USERDRIVINGFORCE_H
#define USERDRIVINGFORCE_H

#include "Tools/Includes.h"
namespace opensim
{

class Settings;
class DrivingForce;
class Temperature;
class PhaseField;
class Node;

class UserDrivingForce : public OPObject
{
 public:

    std::string thisclassname;

    ~UserDrivingForce(void);
    void   Initialize(Settings& locSettings);
    void   SetSpecificDrivingForce(PhaseField& Phase,
            DrivingForce& dGab, int indexA, int indexB, double dfval);
	void   SetSpecificDrivingForceThPhase(PhaseField& Phase,
		DrivingForce& dGab, int indexA, int indexB, double dfval);
    void   CalculateDrivingForce(PhaseField& Phase, Temperature& Tx,
                           DrivingForce& dGab);
    double Energy(PhaseField& Phase, Temperature& Tx);
    double PointEnergy(PhaseField& Phase, Temperature& Tx, int i, int j, int k);
    std::vector<double> LatentHeat;
    int Nphases;
 protected:
 private:

};

}// namespace opensim
#endif
