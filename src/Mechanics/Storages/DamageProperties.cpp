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

#include "Mechanics/Storages/DamageProperties.h"
#include "VTK.h"
#include "PhaseField.h"
#include "Tools.h"
#include "Info.h"
#include "Tools/UserInterface.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "Mechanics/Storages/PlasticProperties.h"
#include "Mechanics/DamageModels/DamageModel.h"
#include "BoundaryConditions.h"
#include "Velocities.h"

namespace opensim
{
using namespace std;


DamageProperties::~DamageProperties(void)
{
    Info::WriteStandard(thisclassname, "Exited normally");
}

void DamageProperties::Initialize(Settings& locSettings)
{
    thisclassname = "DamageProperties";

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    Nphases = locSettings.Nphases;

    Damageflag.Allocate(Nphases);

    EffectiveDamage.Allocate(Nx, Ny, Nz, 1);
    EffectiveDamage_b.Allocate(Nx, Ny, Nz, 1);
    EffectiveDamage_d.Allocate(Nx, Ny, Nz, 1);
    EffectiveDamageN.Allocate(Nx, Ny, Nz, 1);
    peeq.Allocate(Nx, Ny, Nz, 1);

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void DamageProperties::Remesh(const long int newNx, const long int newNy,
                              const long int newNz, BoundaryConditions& BC)
{
    SetBoundaryConditions(BC);
    peeq.Reallocate(newNx, newNy, newNz);

    EffectiveDamageN.Remesh(newNx, newNy, newNz);
    EffectiveDamage.Remesh(newNx, newNy, newNz);
    EffectiveDamage_b.Reallocate(newNx, newNy, newNz);
    EffectiveDamage_d.Reallocate(newNx, newNy, newNz);

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    SetBoundaryConditions(BC);
    Info::WriteStandard(thisclassname, "Remeshed");
}

void DamageProperties::Advect(Velocities& Vel, BoundaryConditions& BC, double dt, int scheme)
{
    if(StrainEquivalentDot.IsNotAllocated())
    {
        StrainEquivalentDot.Allocate(Nx, Ny, Nz, 1);
    }
    if(not StrainEquivalentDot.IsSize(Nx, Ny, Nz))
    {
        StrainEquivalentDot.Reallocate(Nx, Ny, Nz);
    }
    SetBoundaryConditions(BC);

    switch(scheme)
    {
    case Upwind:
    {
        const double dx = Vel.dx;
        const double dx2 = 0.5/dx;
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,StrainEquivalentDot,0,)
        {
            StrainEquivalentDot(i,j,k) = dx2 *
                 ((fabs(Vel.Average(i-1,j,k)[0]) + Vel.Average(i-1,j,k)[0])*peeq(i-1,j,k) +
                  (fabs(Vel.Average(i,j-1,k)[1]) + Vel.Average(i,j-1,k)[1])*peeq(i,j-1,k) +
                  (fabs(Vel.Average(i,j,k-1)[2]) + Vel.Average(i,j,k-1)[2])*peeq(i,j,k-1) +
                  (fabs(Vel.Average(i+1,j,k)[0]) - Vel.Average(i+1,j,k)[0])*peeq(i+1,j,k) +
                  (fabs(Vel.Average(i,j+1,k)[1]) - Vel.Average(i,j+1,k)[1])*peeq(i,j+1,k) +
                  (fabs(Vel.Average(i,j,k+1)[2]) - Vel.Average(i,j,k+1)[2])*peeq(i,j,k+1)) -
                  (fabs(Vel.Average(i,j,k)[0]) +
                   fabs(Vel.Average(i,j,k)[1]) +
                   fabs(Vel.Average(i,j,k)[2])) * peeq(i, j, k)/dx;
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        break;
    }
    case LaxWendroff:
    {
        Info::WriteExit("LaxWendroff advection scheme not supported for DamageStrain",
                        thisclassname, "AdvectEEQ(DP, Vel, BC, dt)");
        exit(13);
        break;
    }
    default:
    {
        Info::WriteExit("Wrong/None advection scheme given in the input file",
                         thisclassname, "AdvectEEQ(DP, Vel, BC, dt)");
        exit(13);
    }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,StrainEquivalentDot,0,)
    {
        for (int n = 0; n < 6; ++n)
        {
            peeq(i, j, k) += StrainEquivalentDot(i, j, k)*dt;
            StrainEquivalentDot(i, j, k) = 0.0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    switch(scheme)
    {
    case Upwind:
    {
        const double dx = Vel.dx;
        const double dx2 = 0.5/dx;
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,StrainEquivalentDot,0,)
        {
            StrainEquivalentDot(i,j,k) = dx2 *
                 ((fabs(Vel.Average(i-1,j,k)[0]) + Vel.Average(i-1,j,k)[0])*EffectiveDamageN(i-1,j,k) +
                  (fabs(Vel.Average(i,j-1,k)[1]) + Vel.Average(i,j-1,k)[1])*EffectiveDamageN(i,j-1,k) +
                  (fabs(Vel.Average(i,j,k-1)[2]) + Vel.Average(i,j,k-1)[2])*EffectiveDamageN(i,j,k-1) +
                  (fabs(Vel.Average(i+1,j,k)[0]) - Vel.Average(i+1,j,k)[0])*EffectiveDamageN(i+1,j,k) +
                  (fabs(Vel.Average(i,j+1,k)[1]) - Vel.Average(i,j+1,k)[1])*EffectiveDamageN(i,j+1,k) +
                  (fabs(Vel.Average(i,j,k+1)[2]) - Vel.Average(i,j,k+1)[2])*EffectiveDamageN(i,j,k+1)) -
                  (fabs(Vel.Average(i,j,k)[0]) +
                   fabs(Vel.Average(i,j,k)[1]) +
                   fabs(Vel.Average(i,j,k)[2])) * EffectiveDamageN(i, j, k)/dx;
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        break;
    }
    case LaxWendroff:
    {
        Info::WriteExit("LaxWendroff advection scheme not supported for DamageStrain",
                        thisclassname, "AdvectEEQ(DP, Vel, BC, dt)");
        exit(13);
        break;
    }
    default:
    {
        Info::WriteExit("Wrong/None advection scheme given in the input file",
                         thisclassname, "AdvectEEQ(DP, Vel, BC, dt)");
        exit(13);
    }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,StrainEquivalentDot,0,)
    {
        for (int n = 0; n < 6; ++n)
        {
            EffectiveDamageN(i, j, k) += StrainEquivalentDot(i, j, k)*dt;
            StrainEquivalentDot(i, j, k) = 0.0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    SetBoundaryConditions(BC);
}

void DamageProperties::SetBoundaryConditions(BoundaryConditions& BC)
{
    BC.SetX(peeq);
    BC.SetY(peeq);
    BC.SetZ(peeq);

    BC.SetX(EffectiveDamageN);
    BC.SetY(EffectiveDamageN);
    BC.SetZ(EffectiveDamageN);

    BC.SetX(EffectiveDamage);
    BC.SetY(EffectiveDamage);
    BC.SetZ(EffectiveDamage);
}

void DamageProperties::SaveDamageN()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EffectiveDamage,0,)
    {
        EffectiveDamageN(i,j,k) = EffectiveDamage(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DamageProperties::CompareDamage()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EffectiveDamage,0,)
    {
        if (EffectiveDamage(i,j,k) < EffectiveDamageN(i,j,k))
        {
            EffectiveDamage(i,j,k) = EffectiveDamageN(i,j,k);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DamageProperties::WriteEffectiveDamageVTK(int tStep)
{
    stringstream buffer;
    std::vector<int> DataTypes {PDScalars};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        VTK::WriteScalar(buffer, EffectiveDamage, "EffectiveDamage");
        VTK::WriteScalar(buffer, peeq, "PEEQ");
        VTK::WriteScalar(buffer, EffectiveDamage_b, "D_b");
        VTK::WriteScalar(buffer, EffectiveDamage_d, "D_d");
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "EffectiveDamage", tStep);
}

std::vector<double> DamageProperties::getMaxDamage()
{
    std::vector<double> result (4, 0);
    for(int i = 0; i < Nx; i++)
    for(int j = 0; j < Ny; j++)
    for(int k = 0; k < Nz; k++)
    if(EffectiveDamage(i,j,k) > result[0])
    {
        result[0] = EffectiveDamage(i,j,k);
        result[1] = i;
        result[2] = j;
        result[3] = k;
    }
    return result;
}

} // namespace openphase
