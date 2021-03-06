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
#include "Mechanics/ElasticityModels/ElasticityKhachaturyan.h"
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

void ElasticityKhachaturyan::CalculateDrivingForce(PhaseField& Phase,
        ElasticProperties& EP, DrivingForce& dGab)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Strains,0,)
    if(Phase.Interface(i,j,k))
    for(auto alpha = Phase.Fields(i, j, k).cbegin();
             alpha < Phase.Fields(i, j, k).cend() - 1; ++alpha)
    for(auto  beta = alpha + 1;
              beta < Phase.Fields(i, j, k).cend();  ++beta)
    {
        double dG_AB = 0.0;		
		
        dMatrix6x6 locCompliance = EP.EffectiveElasticConstants(i,j,k).inverted();
        vStrain ElasticStrains = locCompliance*EP.Stresses(i, j, k);
		
		vStress elasticStresses;
		/*if (EP.NeuberCorrection and (Phase.FieldsStatistics[alpha->index].Phase != Phase.FieldsStatistics[beta->index].Phase))
		{
			elasticStresses = EP.Stresses(i, j, k)*0.01;
		}*/

		if (EP.NeuberCorrection and (Phase.FieldsStatistics[alpha->index].Phase != Phase.FieldsStatistics[beta->index].Phase))
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
			PlasticFlowNeuberMethods::getNeuberDataRO(deviatoricStresses, deviatoricStrains, elasticStresses, dummyStrain);			
			for(int n = 0; n < 3; n++)
		    {
		        elasticStresses[n] += 0.3*StressTrace;
		    }
		}
		else
		{
			elasticStresses = EP.Stresses(i, j, k);
		}
		
        ElasticStrains = locCompliance*elasticStresses;
        
        for(int ii = 0; ii < 6; ii++)
        for(int jj = 0; jj < 6; jj++)
        {
            dG_AB += 0.5*(ElasticStrains[ii])*                                  //elastic strain eps_i
                         (EP.ElasticConstants[beta->index](ii,jj) -
                          EP.ElasticConstants[alpha->index](ii,jj))*            //elastic constants difference
                         (ElasticStrains[jj])                                   //elastic strain eps_j

                       - (ElasticStrains[ii])*                                  //elastic strain eps_i
                         (EP.EffectiveElasticConstants(i,j,k)(ii,jj))*          //effective elastic constant
                         (EP.EigenStrains[beta->index][jj] -
                          EP.EigenStrains[alpha->index][jj]);                   //eigenstrain difference

            /*if(Phase.FieldsStatistics[alpha->index].Stage)
            {
            	dG_AB -= 0.5*(1.0 - Phase.FieldsStatistics[alpha->index].Volume/Phase.RefVolume)*
            			 (EP.EigenStrains[alpha->index][ii] - ElasticStrains[ii])*
						 EP.EffectiveElasticConstants(i,j,k)(ii,jj)*
						 (EP.EigenStrains[alpha->index][jj] - ElasticStrains[jj])*0.0;
            }
            if(Phase.FieldsStatistics[beta->index].Stage)
			{
				dG_AB += 0.5*(1.0 - Phase.FieldsStatistics[beta->index].Volume/Phase.RefVolume)*
						 (EP.EigenStrains[beta->index][ii] - ElasticStrains[ii])*
						 EP.EffectiveElasticConstants(i,j,k)(ii,jj)*
						 (EP.EigenStrains[beta->index][jj] - ElasticStrains[jj])*0.0;
			}*/
        }
        dGab.Raw(i,j,k).add_asym(alpha->index, beta->index, dG_AB);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticityKhachaturyan::CalculateChemicalPotentialContribution(PhaseField& Phase,
                                                          ElasticProperties& EP,
                                              EquilibriumPartitionDiffusion& DF)
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

/*void ElasticityKhachaturyan::CalculateChemicalPotentialContribution(PhaseField& Phase,
                             ElasticProperties& EP, ThermodynamicProperties& TP)
{
    const int Ncomp = EP.Ncomp;
    #pragma omp parallel for collapse(OMP_COLLAPSE_LOOPS) schedule(dynamic, OMP_DYNAMIC_CHUNKSIZE)
    OP_STORAGE_LOOP_ENTIRE(i,j,k,EP.Strains)
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
                                      EP.EffectiveEigenStrains(i,j,k)[ii] +
                                      EP.AppliedStrain[ii]) *
                                     (EP.Strains(i,j,k)[jj] -
                                      EP.EffectiveEigenStrains(i,j,k)[jj] +
                                      EP.AppliedStrain[jj]) +

                                      EP.Lambda({index, comp})[ii]*alpha->value*
                                      EP.EffectiveElasticConstants(i,j,k)(ii,jj)*
                                      (EP.Strains(i,j,k)[jj] -
                                       EP.EffectiveEigenStrains(i,j,k)[jj] +
                                       EP.AppliedStrain[jj]);
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
                locdMuEl += 0.5 * EP.Kappa({index, comp})(ii,jj)*
                                 (EP.Strains(i,j,k)[ii] -
                                  EP.EffectiveEigenStrains(i,j,k)[ii] +
                                  EP.AppliedStrain[ii]) *
                                 (EP.Strains(i,j,k)[jj] -
                                  EP.EffectiveEigenStrains(i,j,k)[jj] +
                                  EP.AppliedStrain[jj]) +

                                  EP.Lambda({index, comp})[ii]*
                                  EP.EffectiveElasticConstants(i,j,k)(ii,jj)*
                                  (EP.Strains(i,j,k)[jj] -
                                   EP.EffectiveEigenStrains(i,j,k)[jj] +
                                   EP.AppliedStrain[jj]);
            }
            TP.ChemicalPotential(i,j,k)({pIndex, comp}) += locdMuEl;
        }
    }
}*/

void ElasticityKhachaturyan::SetEffectiveEigenStrains(PhaseField& Phase,
                                         ElasticProperties& EP, Composition& Cx)
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

void ElasticityKhachaturyan::SetEffectiveEigenStrains(PhaseField& Phase,
                                                      ElasticProperties& EP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveEigenStrains, EP.EffectiveEigenStrains.Bcells(),)
    {
        EP.EffectiveEigenStrains(i,j,k).set_to_zero();
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                EP.EffectiveEigenStrains(i,j,k) +=
                                 EP.EigenStrains[alpha->index]*alpha->value;
            }
        }
        else
        {
            EP.EffectiveEigenStrains(i, j, k) = EP.EigenStrains[Phase.Fields(i,j,k).front().index];
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticityKhachaturyan::SetEffectiveEigenStrains(PhaseField& Phase,
                                        ElasticProperties& EP, Orientations& OR)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveEigenStrains, 0,)
    {
        EP.EffectiveEigenStrains(i,j,k).set_to_zero();
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                EP.EffectiveEigenStrains(i,j,k) +=
                                 EP.EigenStrains[alpha->index]*alpha->value;
            }
        }
        else
        {
            EP.EffectiveEigenStrains(i, j, k) = EP.EigenStrains[Phase.Fields(i,j,k).front().index];
        }
        EP.EffectiveEigenStrains(i,j,k).rotate(OR.Quaternions(i,j,k).RotationMatrix);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticityKhachaturyan::SetEffectiveElasticConstants(PhaseField& Phase,
                                                          ElasticProperties& EP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveElasticConstants, EP.EffectiveElasticConstants.Bcells(),)
    {
        if(Phase.Interface(i,j,k))
        {
            EP.EffectiveElasticConstants(i,j,k).set_to_zero();
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                EP.EffectiveElasticConstants(i,j,k) +=
                             EP.ElasticConstants[alpha->index]*alpha->value;
            }
        }
        else
        {
            EP.EffectiveElasticConstants(i,j,k) =
                               EP.ElasticConstants[Phase.Fields(i,j,k).front().index];
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    EP.AverageElasticConstants.set_to_zero();

    for(unsigned int n = 0; n < Phase.FieldsStatistics.size();n++)
    {
        EP.AverageElasticConstants += EP.ElasticConstants[n]*Phase.FieldsStatistics[n].Volume;
    }
    EP.AverageElasticConstants /= EP.Nx*EP.Ny*EP.Nz;
    EP.AverageCompliences = EP.AverageElasticConstants.inverted();
}

void ElasticityKhachaturyan::SetEffectiveElasticConstants(PhaseField& Phase,
                                        ElasticProperties& EP, Orientations& OR)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveElasticConstants,  EP.EffectiveElasticConstants.Bcells(),)
    {
        if(Phase.Interface(i,j,k))
        {
            EP.EffectiveElasticConstants(i,j,k).set_to_zero();
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                EP.EffectiveElasticConstants(i,j,k) +=
                             EP.ElasticConstants[alpha->index].rotated(
                             OR.Quaternions(i,j,k).RotationMatrix)*alpha->value;
            }
        }
        else
        {
                EP.EffectiveElasticConstants(i,j,k) =
                       EP.ElasticConstants[Phase.Fields(i,j,k).front().index].rotated(
                                          OR.Quaternions(i,j,k).RotationMatrix);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    EP.AverageElasticConstants.set_to_zero();

    for(unsigned int n = 0; n < Phase.FieldsStatistics.size();n++)
    {
        EP.AverageElasticConstants += EP.ElasticConstants[n]*Phase.FieldsStatistics[n].Volume;
    }
    EP.AverageElasticConstants /= EP.Nx*EP.Ny*EP.Nz;
    EP.AverageCompliences = EP.AverageElasticConstants.inverted();
}

void ElasticityKhachaturyan::SetEffectiveElasticConstants(PhaseField& Phase,
                                        ElasticProperties& EP, Orientations& OR,
                                        DamageProperties& DP)
{
    EP.AverageElasticConstants.set_to_zero();
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EP.EffectiveElasticConstants, EP.EffectiveElasticConstants.Bcells(),)
    {
        if(Phase.Interface(i,j,k))
        {
            EP.EffectiveElasticConstants(i,j,k).set_to_zero();
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                EP.EffectiveElasticConstants(i,j,k) +=
                             EP.ElasticConstants[alpha->index].rotated(
                             OR.Quaternions(i,j,k).RotationMatrix)*alpha->value
                             * (1.0 - DP.EffectiveDamage(i,j,k));
            }
        }
        else
        {
                EP.EffectiveElasticConstants(i,j,k) =
                       EP.ElasticConstants[Phase.Fields(i,j,k).front().index].rotated(
                                          OR.Quaternions(i,j,k).RotationMatrix)
                                          * (1.0 - DP.EffectiveDamage(i,j,k));
        }
#ifdef OMP
#pragma omp critical
#endif
        {
            if (i >= 0 and i < EP.Nx and j >= 0 and j < EP.Ny and k >= 0 and k < EP.Nz)
            {
                EP.AverageElasticConstants += EP.EffectiveElasticConstants(i,j,k);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    EP.AverageElasticConstants /= EP.Nx*EP.Ny*EP.Nz;
    EP.AverageCompliences = EP.AverageElasticConstants.inverted();
}

void ElasticityKhachaturyan::SetEffectiveElasticConstants(PhaseField& Phase,
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
    /// TODO: add effect of phase composition to average elastic constants
    for(unsigned int n = 0; n < Phase.FieldsStatistics.size();n++)
    {
        EP.AverageElasticConstants += EP.ElasticConstants[n]*Phase.FieldsStatistics[n].Volume;
    }
    EP.AverageElasticConstants /= EP.Nx*EP.Ny*EP.Nz;
    EP.AverageCompliences = EP.AverageElasticConstants.inverted();
}

}// namespace openphase
