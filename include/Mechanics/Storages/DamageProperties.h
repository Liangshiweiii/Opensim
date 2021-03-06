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

#ifndef DAMAGEPROPERTIES_H
#define DAMAGEPROPERTIES_H

#include "Tools/Includes.h"
namespace opensim
{

class Settings;
class PhaseField;
class Advection;
class ElasticProperties;
class DamageModel;
class Velocities;
class BoundaryConditions;

class DamageProperties : public OPObject                                                       /// Module which stores and handles damage properties
{
 public:
    ~DamageProperties(void);                                                                   /// Destructor, just to indicate that module did not crash at exit

    void Initialize(Settings& locSettings);  /// Initializes the module, allocate the storage, assign internal variables
    void Remesh(const long int newNx, const long int newNy,
                                  const long int newNz, BoundaryConditions& BC);
    void Advect(Velocities& Vel, BoundaryConditions& BC, double dt, int scheme = Upwind);
    void SetBoundaryConditions(BoundaryConditions& BC);
    void SaveDamageN();                                                         /// Saves damage from previous increment
    void CompareDamage();                                                       /// ensures irreversibility of the damage process
    void WriteDamageVTK(int tStep, int Phase);                                  /// Writes strains in VTK format
    void WriteEffectiveDamageVTK(int tStep);                                    /// Writes strains in VTK format
    std::vector<double> getMaxDamage();

    Storage<double>           Damageflag;                                       /// Damage switch for each phase field
    Storage3D<double, 0>      EffectiveDamage;                                  /// Storage for total damage
    Storage3D<double, 0>      EffectiveDamage_b;                                /// Storage for brittle damage
    Storage3D<double, 0>      EffectiveDamage_d;                                /// Storage for ductile damage
    Storage3D<double, 0>      peeq;                                             /// Storage for plastic strain equivalent
    Storage3D<double, 0>      EffectiveDamageN;                                 /// Storage for EffectiveDamage from the previous step


    Storage3D<double, 0>      StrainEquivalentDot;                              /// Storage for advection

    int Nx;
    int Ny;
    int Nz;
    int Nphases;

    bool initialized = false;

 protected:
 private:

};
}// namespace openphase
#endif
