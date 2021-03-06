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

#include "Tools/Node.h"
#include "Tools/NodeMatrix3x3.h"
#include "Tools/Tensor.h"
#include "Info.h"
#include "InterfaceEnergy.h"
#include "InterfaceStress.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "PhaseField.h"

namespace opensim
{

using namespace std;
/******************************************************************************/
InterfaceStress::InterfaceStress(const Settings& locSettings)
{
    this->Initialize(locSettings);
}

void InterfaceStress::Initialize(const Settings& locSettings)
{
    thisclassname = "IsotropicInterfaceStress";
    initialized   = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

const vStress InterfaceStress::operator()(const long int i, const long int j, const long int k,
        const PhaseField& Phase, const InterfaceEnergy& Sigma) const
{
    const double pre1 = 4./Phase.Eta;
    const double pre2 = Phase.Eta * Phase.Eta/Pi;

    vStress locInterfaceStress;
    // Calculate phase gradients and normals
    const NodeV locNormals   = Phase.Normals(i,j,k);
    const NodeV locGradients = Phase.Gradients(i,j,k);

    for (auto it = locNormals.cbegin(); it < locNormals.cend(); ++it)
    {
        // Calculate phase field part of the interface stress
        const double PhiA = abs(Phase.Fields(i,j,k).get(it->indexA));
        const double PhiB = abs(Phase.Fields(i,j,k).get(it->indexB));

        const double interpol = pre1 * (pre2 *
                (locGradients.get(it->indexA) *
                 locGradients.get(it->indexB)) + PhiA * PhiB);

        // Calculate local projection matrix
        dMatrix3x3 Projection;
        Projection.set_to_zero();
        Projection(0,0) = 1 - it->X*it->X;
        Projection(1,0) =   - it->Y*it->X;
        Projection(2,0) =   - it->Z*it->X;
        Projection(0,1) =   - it->X*it->Y;
        Projection(1,1) = 1 - it->Y*it->Y;
        Projection(2,1) =   - it->Z*it->Y;
        Projection(0,2) =   - it->X*it->Z;
        Projection(1,2) =   - it->Y*it->Z;
        Projection(2,2) = 1 - it->Z*it->Z;

        // Calculate local interface stress
        locInterfaceStress -= Projection.VoigtStress() * interpol *
                Sigma(i,j,k,it->indexA,it->indexB);
    }

    return locInterfaceStress;
}

void InterfaceStress::SetEffectiveEigenStrains(ElasticProperties& EP,
        const PhaseField& Phase,
        const InterfaceEnergy& Sigma) const
{

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.EffectiveEigenStrains,0,)
    if (Phase.Interface(i,j,k))
    {
        // Write in effective eigenstrain
        EP.EffectiveEigenStrains(i,j,k) -=
            EP.EffectiveElasticConstants(i,j,k).inverted() *
            (*this)(i,j,k, Phase, Sigma);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

}// namespace opensim
