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

#include "Tools/UserInterface.h"
#include "Compositions.h"
#include "DrivingForce.h"
#include "Diffusion.h"
#include "Info.h"
#include "Mechanics/ElasticityModels/ElasticitySteinbach.h"
#include "Mechanics/PlasticFlow/PlasticFlowNeuberMethods.h"
#include "Mechanics/Storages/DamageProperties.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "Orientations.h"
#include "PhaseField.h"
#include "Settings.h"
#include "VTK.h"

namespace opensim
{

using namespace std;

void ElasticitySteinbach::CalculateDrivingForce(PhaseField& Phase,
                                          ElasticProperties& EP,
                                          DrivingForce& dGab)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Stresses, 0,)
    if(Phase.Interface(i,j,k))
    {
        dMatrix6x6 locCompliance = EP.EffectiveElasticConstants(i,j,k).inverted();
        vStrain ElasticStrains = locCompliance*EP.Stresses(i, j, k);
        
		vStress elasticStresses;		
		if (EP.NeuberCorrection)
		{
			vStress deviatoricStresses;
		    vStrain deviatoricStrains;
		    vStrain dummyStrain;
		    double StressTrace = EP.Stresses(i, j, k).trace();
		    double StrainTrace = ElasticStrains.trace();
		    for(int n = 0; n < 3; n++)
		    {
		        deviatoricStresses[n] = EP.Stresses(i, j, k)[n] - 0.3*StressTrace;
		        deviatoricStrains[n] = ElasticStrains[n] - 0.3*StrainTrace;
		    }		    
			PlasticFlowNeuberMethods::getNeuberData(deviatoricStresses, deviatoricStrains, elasticStresses, dummyStrain);			
			for(int n = 0; n < 3; n++)
		    {
		        elasticStresses[n] += 0.3*StressTrace;
		    }
		}
		else
		{
			elasticStresses = EP.Stresses(i, j, k);
		}

        for(auto alpha = Phase.Fields(i, j, k).cbegin();
                 alpha < Phase.Fields(i, j, k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta < Phase.Fields(i, j, k).cend();  ++beta)
        {
            double dG_AB = 0.0;

            for(int ii = 0; ii < 6; ii++)
            {
                dG_AB += elasticStresses[ii]*
                        (EP.EigenStrains[alpha->index][ii] -
                         EP.EigenStrains[beta->index][ii]);

                for(int jj = 0; jj < 6; jj++)
                {
                    dG_AB += 0.5*elasticStresses[ii]*
                                (EP.Compliences[alpha->index](ii,jj) -
                                 EP.Compliences[beta->index](ii,jj))*
								elasticStresses[jj];
                }
            }
            dGab.Raw(i,j,k).add_asym(alpha->index, beta->index, dG_AB);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

/*void ElasticitySteinbach::CalculateChemicalPotentialContribution(PhaseField& Phase,
                           ElasticProperties& EP, ThermodynamicFunctionsOLD& TF)
{
    const int Nphases = EP.Nphases;
    const int Ncomp = EP.Ncomp;
    #pragma omp parallel for collapse(OMP_COLLAPSE_LOOPS) schedule(dynamic, OMP_DYNAMIC_CHUNKSIZE)
    //OP_STORAGE_LOOP_ENTIRE(i,j,k,EP.Strains)
    OP_STORAGE_LOOP_INTERIOR(i,j,k,EP.Strains)
    for(int comp = 0; comp < Ncomp; comp++)
    {
        for (int n = 0; n < Nphases; n++)
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                double locdMuEl = 0.0;
                int index = alpha->index;
                int pIndex = Phase.FieldsStatistics[index].Phase;
                for(int ii = 0; ii < 6; ii++)
                {
                    locdMuEl -= EP.Lambda({index, comp})[ii]*alpha->value*
                                EP.Stresses(i, j, k)[ii];
                    for(int jj = 0; jj < 6; jj++)
                    {
                        locdMuEl += 0.5 * EP.Kappa({index, comp})(ii,jj)*alpha->value*
                                         (EP.Strains(i,j,k)[ii] -
                                          EP.EffectiveEigenStrains(i,j,k)[ii]) *
                                         (EP.Strains(i,j,k)[jj] -
                                          EP.EffectiveEigenStrains(i,j,k)[jj]);
                    }
                }
                TF.ChemicalPotential(i,j,k)({pIndex,comp}) += locdMuEl;
            }
        }
        else
        {
            double locdMuEl = 0.0;
            int index = Phase.Fields(i,j,k).front().index;
            int pIndex = Phase.FieldsStatistics[index].Phase;
            for(int ii = 0; ii < 6; ii++)
            {
                locdMuEl -= EP.Lambda({index, comp})[ii]*EP.Stresses(i, j, k)[ii];
                for(int jj = 0; jj < 6; jj++)
                {
                    locdMuEl += 0.5 * EP.Kappa({index, comp})(ii,jj)*
                                     (EP.Strains(i,j,k)[ii] -
                                      EP.EffectiveEigenStrains(i,j,k)[ii]) *
                                     (EP.Strains(i,j,k)[jj] -
                                      EP.EffectiveEigenStrains(i,j,k)[jj]);
                }
            }
            TF.ChemicalPotential(i,j,k)({pIndex,comp}) += locdMuEl;
        }
    }
}*/

void ElasticitySteinbach::CalculateChemicalPotentialContribution(PhaseField& Phase,
                                                        ElasticProperties& EP,
                                            Diffusion& DF)
{
    const int Nphases = EP.Nphases;
    const int Ncomp = EP.Ncomp;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Strains, EP.Strains.Bcells(),)
    for(int comp = 0; comp < Ncomp; comp++)
    {
        for (int n = 0; n < Nphases; n++) DF.dMu(i,j,k)({n}) = 0.0;
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                double locdMuEl = 0.0;
                int index = alpha->index;
                int pIndex = Phase.FieldsStatistics[index].Phase;
                for(int ii = 0; ii < 6; ii++)
                {
                    locdMuEl -= EP.Lambda({index, comp})[ii]*alpha->value*
                                EP.Stresses(i, j, k)[ii];
                    for(int jj = 0; jj < 6; jj++)
                    {
                        locdMuEl += 0.5 * EP.Kappa({index, comp})(ii,jj)*alpha->value*
                                         (EP.Strains(i,j,k)[ii] -
                                          EP.EffectiveEigenStrains(i,j,k)[ii]) *
                                          EP.ElasticConstants[index](ii,jj) *
                                         (EP.Strains(i,j,k)[jj] -
                                          EP.EffectiveEigenStrains(i,j,k)[jj]);
                    }
                }
                DF.dMu(i,j,k)({pIndex}) += locdMuEl;
            }
        }
        else
        {
            double locdMuEl = 0.0;
            int index = Phase.Fields(i,j,k).front().index;
            int pIndex = Phase.FieldsStatistics[index].Phase;
            for(int ii = 0; ii < 6; ii++)
            {
                locdMuEl -= EP.Lambda({index, comp})[ii]*EP.Stresses(i, j, k)[ii];
                for(int jj = 0; jj < 6; jj++)
                {
                    locdMuEl += 0.5 * EP.Kappa({index, comp})(ii,jj)*
                                     (EP.Strains(i,j,k)[ii] -
                                      EP.EffectiveEigenStrains(i,j,k)[ii]) *
                                      EP.ElasticConstants[index](ii,jj) *
                                     (EP.Strains(i,j,k)[jj] -
                                      EP.EffectiveEigenStrains(i,j,k)[jj]);
                }
            }
            DF.dMu(i,j,k)({pIndex}) = locdMuEl;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

/*void ElasticitySteinbach::CalculateChemicalPotentialContrib(PhaseField& Phase,
                            ElasticProperties& EP, ThermodynamicProperties& TP)
{
    const int Ncomp = EP.Ncomp;
    #pragma omp parallel for collapse(OMP_COLLAPSE_LOOPS) schedule(dynamic, OMP_DYNAMIC_CHUNKSIZE)
    //OP_STORAGE_LOOP_ENTIRE(i,j,k,EP.Strains)
    OP_STORAGE_LOOP_INTERIOR(i,j,k,EP.Strains)
    for(int comp = 0; comp < Ncomp; comp++)
    {
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                double locdMuEl = 0.0;
                int index = alpha->index;
                int pIndex = Phase.FieldsStatistics[index].Phase;

                for(int ii = 0; ii < 6; ii++)
                for(int jj = 0; jj < 6; jj++)
                {
                    locdMuEl += 0.5 * EP.Kappa({index, comp})(ii,jj)*alpha->value*
                                     (EP.Strains(i,j,k)[ii] -
                                      EP.EffectiveEigenStrains(i,j,k)[ii]) *
                                     (EP.Strains(i,j,k)[jj] -
                                      EP.EffectiveEigenStrains(i,j,k)[jj]) +

                                      EP.Lambda({index, comp})[ii]*alpha->value*
                                      EP.EffectiveElasticConstants(i,j,k)(ii,jj)*
                                      (EP.Strains(i,j,k)[jj] -
                                       EP.EffectiveEigenStrains(i,j,k)[jj]);
                }
                TP.ChemicalPotential(i,j,k)({pIndex, comp}) += locdMuEl;
            }
        }
        else
        {
            double locdMuEl = 0.0;
            int index = Phase.Fields(i,j,k).front().index;
            int pIndex = Phase.FieldsStatistics[index].Phase;

            for(int ii = 0; ii < 6; ii++)
            for(int jj = 0; jj < 6; jj++)
            {
                locdMuEl += 0.5 * EP.Kappa({index,comp})(ii,jj)*
                                 (EP.Strains(i,j,k)[ii] -
                                  EP.EffectiveEigenStrains(i,j,k)[ii]) *
                                 (EP.Strains(i,j,k)[jj] -
                                  EP.EffectiveEigenStrains(i,j,k)[jj]) +

                                  EP.Lambda({index, comp})[ii]*
                                  EP.EffectiveElasticConstants(i,j,k)(ii,jj)*
                                  (EP.Strains(i,j,k)[jj] -
                                   EP.EffectiveEigenStrains(i,j,k)[jj]);
            }
            TP.ChemicalPotential(i,j,k)({pIndex, comp}) += locdMuEl;
        }
    }
}*/

void ElasticitySteinbach::SetEffectiveEigenStrains(PhaseField& Phase,
                                                   ElasticProperties& EP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveEigenStrains, EP.EffectiveEigenStrains.Bcells(),)
    {
        EP.EffectiveEigenStrains(i,j,k).set_to_zero();
        if(!(Phase.Interface(i,j,k)))
        {
            EP.EffectiveEigenStrains(i,j,k) =
                                 EP.EigenStrains[Phase.Fields(i,j,k).front().index];
        }
        else
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                EP.EffectiveEigenStrains(i,j,k) +=
                                     EP.EigenStrains[alpha->index]*alpha->value;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySteinbach::SetEffectiveEigenStrains(PhaseField& Phase,
                                        ElasticProperties& EP, Orientations& OR)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveEigenStrains, EP.EffectiveEigenStrains.Bcells(),)
    {
        EP.EffectiveEigenStrains(i,j,k).set_to_zero();
        if(!(Phase.Interface(i,j,k)))
        {
            EP.EffectiveEigenStrains(i,j,k) = EP.EigenStrains[Phase.Fields(i,j,k).front().index];
        }
        else
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                EP.EffectiveEigenStrains(i, j, k) +=
                                     EP.EigenStrains[alpha->index]*alpha->value;
            }
        }
        EP.EffectiveEigenStrains(i,j,k).rotate(OR.Quaternions(i,j,k).RotationMatrix);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySteinbach::SetEffectiveEigenStrains(PhaseField& Phase, ElasticProperties& EP, Composition& Cx)
{
    const int Ncomp = EP.Ncomp;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveEigenStrains, EP.EffectiveEigenStrains.Bcells(),)
    {
        EP.EffectiveEigenStrains(i,j,k).set_to_zero();
        if(!(Phase.Interface(i,j,k)))
        {
            int index = Phase.Fields(i,j,k).front().index;
            int pIndex = Phase.FieldsStatistics[index].Phase;
            EP.EffectiveEigenStrains(i,j,k) = EP.EigenStrains[index];

            for(int comp = 0; comp < Ncomp; comp++)
            {
                double delta = (Cx.Phase(i,j,k)({pIndex, comp}) -
                                                   EP.Cref({pIndex, comp}));
                for(int dir = 0; dir < 6; dir++)
                {
                    EP.EffectiveEigenStrains(i,j,k)[dir] += delta *
                                              EP.Lambda({index, comp})[dir];
                }
            }
        }
        else
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                int index = alpha->index;
                int pIndex = Phase.FieldsStatistics[index].Phase;
                double phiAlpha = alpha->value;

                EP.EffectiveEigenStrains(i,j,k) +=
                                            EP.EigenStrains[index]*phiAlpha;

                for(int comp = 0; comp < Ncomp; comp++)
                {
                    double delta = (Cx.Phase(i,j,k)({pIndex, comp}) -
                                                   EP.Cref({pIndex, comp}));

                    for(int dir = 0; dir < 6; dir++)
                    {
                        EP.EffectiveEigenStrains(i,j,k)[dir] += delta*
                                     phiAlpha*EP.Lambda({index, comp})[dir];
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySteinbach::SetEffectiveEigenStrains(PhaseField& Phase,
                       ElasticProperties& EP, Composition& Cx, Orientations& OR)
{
    const int Ncomp = EP.Ncomp;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveEigenStrains, EP.EffectiveElasticConstants.Bcells(),)
    {
        EP.EffectiveEigenStrains(i,j,k).set_to_zero();
        if(!(Phase.Interface(i,j,k)))
        {
            int index = Phase.Fields(i,j,k).front().index;
            int pIndex = Phase.FieldsStatistics[index].Phase;
            EP.EffectiveEigenStrains(i,j,k) = EP.EigenStrains[index];

            for(int comp = 0; comp < Ncomp; comp++)
            {
                double delta = (Cx.Phase(i,j,k)({pIndex, comp}) -
                                                   EP.Cref({pIndex, comp}));
                for(int dir = 0; dir < 6; dir++)
                {
                    EP.EffectiveEigenStrains(i,j,k)[dir] += delta *
                                              EP.Lambda({index, comp})[dir];
                }
            }
        }
        else
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                int index = alpha->index;
                int pIndex = Phase.FieldsStatistics[index].Phase;
                double phiAlpha = alpha->value;

                EP.EffectiveEigenStrains(i,j,k) +=
                                            EP.EigenStrains[index]*phiAlpha;

                for(int comp = 0; comp < Ncomp; comp++)
                {
                    double delta = (Cx.Phase(i,j,k)({pIndex, comp}) -
                                                   EP.Cref({pIndex, comp}));

                    for(int dir = 0; dir < 6; dir++)
                    {
                        EP.EffectiveEigenStrains(i,j,k)[dir] += delta*
                                     phiAlpha*EP.Lambda({index, comp})[dir];
                    }
                }
            }
        }
        EP.EffectiveEigenStrains(i,j,k).rotate(OR.Quaternions(i,j,k).RotationMatrix);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySteinbach::SetEffectiveElasticConstants(PhaseField& Phase,
                                                       ElasticProperties& EP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveElasticConstants, EP.EffectiveElasticConstants.Bcells(),)
    {
        EP.EffectiveElasticConstants(i,j,k).set_to_zero();
        if(Phase.Interface(i,j,k))
        {
            dMatrix6x6 TempCompliences;
            TempCompliences.set_to_zero();

            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                TempCompliences += (EP.Compliences[alpha->index]*alpha->value);
            }
            EP.EffectiveElasticConstants(i,j,k) = TempCompliences.inverted();
        }
        else
        {
            EP.EffectiveElasticConstants(i,j,k) =
              EP.ElasticConstants[Phase.Fields(i,j,k).front().index];
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    EP.AverageElasticConstants.set_to_zero();
    EP.AverageCompliences.set_to_zero();

    for(unsigned int n = 0; n < Phase.FieldsStatistics.size();n++)
    {
        EP.AverageElasticConstants += EP.ElasticConstants[n]*Phase.FieldsStatistics[n].Volume;
        EP.AverageCompliences += EP.ElasticConstants[n].inverted()*Phase.FieldsStatistics[n].Volume;
    }
    EP.AverageElasticConstants /= EP.Nx*EP.Ny*EP.Nz;
    EP.AverageCompliences /= EP.Nx*EP.Ny*EP.Nz;
}

void ElasticitySteinbach::SetEffectiveElasticConstants(PhaseField& Phase,
                                        ElasticProperties& EP, Orientations& OR)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveElasticConstants, EP.EffectiveElasticConstants.Bcells(),)
    {
        EP.EffectiveElasticConstants(i,j,k).set_to_zero();
        if(Phase.Interface(i,j,k))
        {
            dMatrix6x6 TempCompliences;
            TempCompliences.set_to_zero();

            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                dMatrix6x6 tempStiff =
                   EP.ElasticConstants[alpha->index].rotated(OR.Quaternions(i,j,k).RotationMatrix);
                TempCompliences += (tempStiff.inverted()*alpha->value);
            }
            EP.EffectiveElasticConstants(i,j,k) = TempCompliences.inverted();
        }
        else
        {
            EP.EffectiveElasticConstants(i,j,k) =
               EP.ElasticConstants[Phase.Fields(i,j,k).front().index].rotated(OR.Quaternions(i,j,k).RotationMatrix);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    EP.AverageElasticConstants.set_to_zero();
    EP.AverageCompliences.set_to_zero();

    for(unsigned int n = 0; n < Phase.FieldsStatistics.size();n++)
    {
        EP.AverageElasticConstants += EP.ElasticConstants[n]*Phase.FieldsStatistics[n].Volume;
        EP.AverageCompliences += EP.ElasticConstants[n].inverted()*Phase.FieldsStatistics[n].Volume;
    }
    EP.AverageElasticConstants /= EP.Nx*EP.Ny*EP.Nz;
    EP.AverageCompliences /= EP.Nx*EP.Ny*EP.Nz;
}

void ElasticitySteinbach::SetEffectiveElasticConstants(PhaseField& Phase,
                                        ElasticProperties& EP, Orientations& OR,
                                        DamageProperties& DP)
{
    EP.AverageElasticConstants.set_to_zero();
    EP.AverageCompliences.set_to_zero();

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveElasticConstants, EP.EffectiveElasticConstants.Bcells(),)
    {
        EP.EffectiveElasticConstants(i,j,k).set_to_zero();
        if(Phase.Interface(i,j,k))
        {
            dMatrix6x6 TempCompliences;
            TempCompliences.set_to_zero();

            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                dMatrix6x6 tempStiff =
                   EP.ElasticConstants[alpha->index].rotated(OR.Quaternions(i,j,k).RotationMatrix)
                                                               * (1.0 - DP.EffectiveDamage(i,j,k));
                TempCompliences += (tempStiff.inverted()*alpha->value);
            }
            EP.EffectiveElasticConstants(i,j,k) = TempCompliences.inverted();
        }
        else
        {
            EP.EffectiveElasticConstants(i,j,k) =
               EP.ElasticConstants[Phase.Fields(i,j,k).front().index].rotated(OR.Quaternions(i,j,k).RotationMatrix)
                                                                       * (1.0 - DP.EffectiveDamage(i,j,k));
        }
#ifdef OMP
#pragma omp critical
#endif
        {
            if (i >= 0 and i < EP.Nx and j >= 0 and j < EP.Ny and k >= 0 and k < EP.Nz)
            {
                EP.AverageElasticConstants += EP.EffectiveElasticConstants(i,j,k);
                EP.AverageCompliences += EP.EffectiveElasticConstants(i,j,k).inverted();
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    EP.AverageElasticConstants /= EP.Nx*EP.Ny*EP.Nz;
    EP.AverageCompliences /= EP.Nx*EP.Ny*EP.Nz;
}

void ElasticitySteinbach::SetEffectiveElasticConstants(PhaseField& Phase,
                                         ElasticProperties& EP, Composition& Cx)
{
    const int Ncomp = EP.Ncomp;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveElasticConstants, EP.EffectiveElasticConstants.Bcells(),)
    {
        EP.EffectiveElasticConstants(i,j,k).set_to_zero();
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                double PhiAlpha = alpha->value;
                int index = alpha->index;
                int pIndex = Phase.FieldsStatistics[index].Phase;

                for(int ii = 0; ii < 6; ii++)
                for(int jj = 0; jj < 6; jj++)
                {
                    double deltaCij = 0.0;
                    for(int comp = 0; comp < Ncomp; comp ++)
                    {
                        deltaCij += EP.Kappa({index, comp})(ii,jj)*
                        (Cx.Phase(i,j,k)({pIndex, comp}) - EP.Cref({pIndex, comp}));
                    }
                    EP.EffectiveElasticConstants(i,j,k)(ii,jj) += PhiAlpha *
                            (EP.ElasticConstants[alpha->index](ii,jj) + deltaCij);
                }
            }
        }
        else
        {
            int index = Phase.Fields(i,j,k).front().index;
            int pIndex = Phase.FieldsStatistics[index].Phase;

            for(int ii = 0; ii < 6; ii++)
            for(int jj = 0; jj < 6; jj++)
            {
                double deltaCij = 0.0;
                for(int comp = 0; comp < Ncomp; comp ++)
                {
                    deltaCij += EP.Kappa({index, comp})(ii,jj)*
                    (Cx.Phase(i,j,k)({pIndex, comp}) - EP.Cref({pIndex, comp}));
                }
                EP.EffectiveElasticConstants(i,j,k)(ii,jj) =
                               EP.ElasticConstants[index](ii,jj) + deltaCij;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    EP.AverageElasticConstants.set_to_zero();
    EP.AverageCompliences.set_to_zero();
    /// TODO: add effect of phase composition to average elastic constants
    for(unsigned int n = 0; n < Phase.FieldsStatistics.size();n++)
    {
        EP.AverageElasticConstants += EP.ElasticConstants[n]*Phase.FieldsStatistics[n].Volume;
        EP.AverageCompliences += EP.ElasticConstants[n].inverted()*Phase.FieldsStatistics[n].Volume;
    }
    EP.AverageElasticConstants /= EP.Nx*EP.Ny*EP.Nz;
    EP.AverageCompliences /= EP.Nx*EP.Ny*EP.Nz;
}

}// namespace openphase
