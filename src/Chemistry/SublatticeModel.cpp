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

#include "Info.h"
#include "Chemistry/SublatticeModel.h"

namespace opensim
{

using namespace std;

void SublatticeModel::Initialize(void)
{
    Ncons = getNcons();
    hasVacancies = gethasVacancies();
}

int SublatticeModel::getNcons()
{
    /**This will return the number of constituents of this sublattice.*/
    return Constituent.size();    
}

bool SublatticeModel::gethasVacancies()
{
    /**This function will loop over all constituents and checks, if one of them
    is considered a vacancy.*/
    bool temp = false;
    for(int i = 0; i < Ncons; i++)
    if(Constituent[i].Index < 0)
    {
        temp = true;
    }
    return temp;
}

bool SublatticeModel::isElementPresent(int i)
{
    /**This function will loop over all constituents and checks, if element i
    is present on this sublattice. "i" has to be the index, not the component
    number!!*/
    bool temp = false;
    for(int n = 0; n < Ncons; n++)
    {
        if(Constituent[n].Index == i)
        temp = true;
    }
    return temp;
}

SublatticeModel::SublatticeModel(double n_sites, unsigned int n_cons)
{
    /**Constructor*/
    Site = n_sites;
}

}
