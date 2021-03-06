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

//#include "UserDrivingForce.h"
#include "Info.h"
#include "Settings.h"
#include "PhaseField.h"
#include "DrivingForce.h"
#include "Tools/Node.h"
#include "Temperatures.h"

namespace opensim
{
using namespace std;

void UserDrivingForce::Initialize(Settings& locSettings)
{
    thisclassname = "UserDrivingForce";

    Nphases = locSettings.Nphases;
    LatentHeat.resize(Nphases);
    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

UserDrivingForce::~UserDrivingForce(void)
{
    Info::WriteStandard(thisclassname, "Exited normally");
}

void UserDrivingForce::SetSpecificDrivingForce(PhaseField& Phase,
        DrivingForce& dGab, int indexA, int indexB, double dfval)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Interface(i,j,k))
		{
			dGab.Raw(i, j, k).clear();
            bool foundindexA = Phase.PhaseFieldPresent(i,j,k,indexA);
            bool foundindexB = Phase.PhaseFieldPresent(i,j,k,indexB);

            if (foundindexA == true and foundindexB == true)
            {
                dGab.Raw(i,j,k).add_asym(indexA, indexB, dfval);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void UserDrivingForce::SetSpecificDrivingForceThPhase(PhaseField& Phase,
                       DrivingForce& dGab, int indexA, int indexB, double dfval)
    {
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields, 0, )
    {
        if (Phase.Interface(i, j, k))
        {
            dGab.Raw(i, j, k).clear();
            for (auto alpha = Phase.Fields(i, j, k).cbegin();
                      alpha != Phase.Fields(i, j, k).cend(); ++alpha)
            for (auto beta = Phase.Fields(i, j, k).cbegin();
                      beta != Phase.Fields(i, j, k).cend(); ++beta)
            {
                int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                int pIndexB = Phase.FieldsStatistics[beta->index].Phase;

                if (indexA == pIndexA and indexB == pIndexB)
                {
                    dGab.Raw(i,j,k).add_asym(alpha->index, beta->index, dfval);
                }
                else if(indexB == pIndexA and indexA == pIndexB)
                {
                    dGab.Raw(i,j,k).add_asym(beta->index, alpha->index, dfval);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void UserDrivingForce::CalculateDrivingForce(PhaseField& Phase,
                                        Temperature& Tx,
                                        DrivingForce& dGab)
{               
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend()-1; ++alpha)
            for(auto beta = alpha + 1;
                     beta < Phase.Fields(i,j,k).cend(); ++beta)
            {
                int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;
                double dG_AB = (Tx(i,j,k)-Tx.T0)/Tx.T0*(LatentHeat[pIndexB] -
                                                        LatentHeat[pIndexA]); 
                dGab.Raw(i,j,k).add_asym(alpha->index, beta->index, dG_AB);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double UserDrivingForce::Energy(PhaseField& Phase, Temperature& Tx)
{
    double Energy = 0.0;

    const int Nx = Phase.Nx;
    const int Ny = Phase.Ny;
    const int Nz = Phase.Nz;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        Energy += PointEnergy(Phase, Tx, i, j, k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    
    return Energy/double(Nx*Ny*Nz);
}

double UserDrivingForce::PointEnergy(PhaseField& Phase, Temperature& Tx,
                                      int i, int j, int k)
{
    double Energy = 0.0;

    for(auto alpha = Phase.Fields(i,j,k).cbegin();
             alpha < Phase.Fields(i,j,k).cend(); ++alpha)
    {
        int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
        Energy += (Tx(i,j,k)-Tx.T0)/Tx.T0*alpha->value*LatentHeat[pIndexA];
    }
    return Energy;
}
}// namespace opensim
