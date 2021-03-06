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

#ifndef PLASTICPROPERTIES_H
#define PLASTICPROPERTIES_H

#include "Tools/Includes.h"
namespace opensim
{

class Settings;
class Composition;
class PhaseField;
class Advection;
class BoundaryConditions;
class Temperature;
class Orientations;
class Velocities;

class PlasticProperties : public OPObject                                   /// Module which stores and handles data storages for the plastic flow rule model
{
public:
    PlasticProperties(){};
    PlasticProperties(Settings& locSettings);
    ~PlasticProperties();

    void Initialize(Settings& locSettings);
    using OPObject::ReadInput;

    void SetBoundaryConditions(const BoundaryConditions& BC);
    void Remesh(int newNx, int newNy, int newNz, BoundaryConditions& BC);
    void Write(int tStep, bool legacy_format = true);
    void Read(const BoundaryConditions& BC, int tStep, bool legacy_format = true);
    void Read(std::string FileName, bool legacy_format);

    void Advect(Velocities& Vel, BoundaryConditions& BC,
                            double dt, int scheme = Upwind);

    void WritePlasticStrainVTK(int tStep);

    int Nx;
    int Ny;
    int Nz;

    Storage3D<vStrain, 0> PlasticStrain;
    Storage3D<vStrain, 0> StrainsDot;                                           /// storage for the increment of the strains due to advection

protected:
private:

};
}// namespace openphase
#endif
