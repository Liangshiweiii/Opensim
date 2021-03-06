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

#include <Mechanics/Storages/PlasticProperties.h>
#include "Mechanics/LargeDeformations/LargeDeformationsMethods.h"
#include "Info.h"
#include "Velocities.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "Orientations.h"
#include "Settings.h"
#include "DrivingForce.h"
#include "PhaseField.h"
#include "Temperatures.h"
#include "Compositions.h"
#include "BoundaryConditions.h"
#include "InterfaceEnergy.h"
#include "InterfaceMobility.h"
#include "Tools/Quaternion.h"
#include "Mechanics/Storages/DamageProperties.h"
#include "Mechanics/PlasticFlow/PlasticFlowCP.h"
#include "VTK.h"

namespace opensim
{
using namespace std;

void LargeDeformationsMethods::CalculateVelocityGradient(ElasticProperties& EP,
                                                                Velocities& Vel)
{
    double dx2Inv = 0.5/EP.dx;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Vel.Average, 0,)
    for(int dir = 0; dir < 3; dir++)
    {
        EP.VelocityGradient(i,j,k)(dir,0) = dx2Inv*
                        (Vel.Average(i+1,j,k)[dir] - Vel.Average(i-1,j,k)[dir]);

        EP.VelocityGradient(i,j,k)(dir,1) = dx2Inv*
                        (Vel.Average(i,j+1,k)[dir] - Vel.Average(i,j-1,k)[dir]);

        EP.VelocityGradient(i,j,k)(dir,2) = dx2Inv*
                        (Vel.Average(i,j,k+1)[dir] - Vel.Average(i,j,k-1)[dir]);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

/*
void LargeDeformations::CalculateAngularVelocity(ElasticProperties& EP)
{
    int Nx = EP.Nx;
    int Ny = EP.Ny;
    int Nz = EP.Nz;

    #pragma omp parallel // OMP BEGIN
    {
        int nThreads = omp_get_num_threads();
        int myThread = omp_get_thread_num();

        int BegR = (myThread*Nx)/nThreads;
        int EndR = ((myThread+1)*Nx)/nThreads;

        for(int i = BegR+1; i < EndR+1; i++)
        for(int j =    1; j <   Ny+1; j++)
        for(int k =    1; k <   Nz+1; k++)
        {
            //dMatrix3x3 L = EP.VelocityGradient(i,j,k);
            //EP.AngularVelocity(i,j,k)=((L-L.transposed())*0.5);

        }
    }// OMP END
}*/

void LargeDeformationsMethods::AccumulateStresses(ElasticProperties& EP,
                                                                      double dt)
{
OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Stresses, 0,)
    {
        dMatrix3x3 SigmaAcc = EP.Stresses(i,j,k).tensor();
        dMatrix3x3 L = EP.VelocityGradient(i,j,k);
        EP.Stresses(i,j,k) += EP.StressIncrements(i,j,k) +
            (SigmaAcc*L.transposed() + L*SigmaAcc - SigmaAcc*L.trace()).VoigtStress()*dt;
        // To avoid second acc. in case of remeshing
        EP.StressIncrements(i,j,k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::ConsiderStiffnessIncrement(ElasticProperties& EP,
                                        Storage3D<dMatrix6x6, 0>& StiffnessOLD)
{
    cout<<"entering consider stiffness"<<endl;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.EffectiveElasticConstants, 0,)
    if(StiffnessOLD(i,j,k)(0,0) > 0.0)
    {
        dMatrix6x6 dC = EP.EffectiveElasticConstants(i,j,k)-StiffnessOLD(i,j,k);
        EP.EffectiveEigenStrains(i,j,k) +=
                         (dC.inverted()*EP.Stresses(i,j,k)).Ln();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    cout<<"finished"<<endl;
}

void LargeDeformationsMethods::AccumulateStressesRotation(
                              ElasticProperties& EP, Orientations& OR,
                              Storage3D<dMatrix3x3, 0>& RotationsOLD,
                          /*Storage3D<dMatrix6x6, 0>& StiffnessOLD,*/ double dt)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Strains, 0,)
    {
        dMatrix3x3 R;
//        if(RotationsOLD(i,j,k).norm() < 1.)
//        {
//            cout<<"at "<<i<<"|"<<j<<"|"<<k<<endl;
//            cout<<"norm = "<<RotationsOLD(i,j,k).norm()<<endl;
//        }
//        if(StiffnessOLD(i,j,k)(0,0)!=0.0)
//        {
//            const dMatrix6x6 dC =
//                    EP.EffectiveElasticConstants(i,j,k)-StiffnessOLD(i,j,k);
//            const vStrain epsilon = StiffnessOLD(i,j,k).inverted()*EP.Stresses(i,j,k);
//            EP.StressIncrements(i,j,k) += dC*epsilon;
//        }
        R = OR.Quaternions(i,j,k).RotationMatrix*
               RotationsOLD(i,j,k).inverted();
        R /= pow(R.determinant(), 1./3.);
#ifdef DEBUG
        if(fabs(R.determinant()-1.) > 1e-8)
        {
            cout<<"Rotation increment at ("<<i<<", "<<j<<", "<<k
                <<") is not stretchfree!\n"
                <<"|dR| - 1 = "<<R.determinant() - 1.<<endl
                <<"|R_OLD| - 1 = "<<RotationsOLD(i,j,k).determinant() - 1.<<endl
                <<"|R_NEW| - 1 = "<<OR.Quaternions(i,j,k).RotationMatrix.determinant() - 1.<<endl;

            cout<<"before correction:"<<endl<<"R:"<<endl<<R.print()<<endl;
            cout<<"correction factor: "<< pow(R.determinant(), 1./3.)<<endl;
            R /= pow(R.determinant(), 1./3.);
            cout<<"after correction"<<endl<<"R:"<<endl<<R.print()<<endl;
            cout<<"|dR|-1 = "<<R.determinant() - 1.<<endl;
//            exit(13);
        }
#endif
        EP.Stresses(i,j,k).rotate(R);
        EP.Stresses(i,j,k) += EP.StressIncrements(i,j,k);
        //To avoid second acc. in case of remeshing
        EP.StressIncrements(i,j,k).set_to_zero();
        EP.Strains(i,j,k) += EP.StrainIncrements(i,j,k).Ln();
		EP.StrainIncrements(i,j,k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}



void LargeDeformationsMethods::RotatePlasticStrains(
                              PlasticProperties& PFP, Orientations& OR,
                              Storage3D<dMatrix3x3, 0>& RotationsOLD,
                          /*Storage3D<dMatrix6x6, 0>& StiffnessOLD,*/ double dt)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, PFP.PlasticStrain, 0,)
    {
        dMatrix3x3 R;
//        if(RotationsOLD(i,j,k).norm() < 1.)
//        {
//            cout<<"at "<<i<<"|"<<j<<"|"<<k<<endl;
//            cout<<"norm = "<<RotationsOLD(i,j,k).norm()<<endl;
//        }
//        if(StiffnessOLD(i,j,k)(0,0)!=0.0)
//        {
//            const dMatrix6x6 dC =
//                    EP.EffectiveElasticConstants(i,j,k)-StiffnessOLD(i,j,k);
//            const vStrain epsilon = StiffnessOLD(i,j,k).inverted()*EP.Stresses(i,j,k);
//            EP.StressIncrements(i,j,k) += dC*epsilon;
//        }
        R = OR.Quaternions(i,j,k).RotationMatrix*
               RotationsOLD(i,j,k).inverted();
        R /= pow(R.determinant(), 1./3.);
#ifdef DEBUG
        if(fabs(R.determinant()-1.) > 1e-8)
        {
            cout<<"Rotation increment at ("<<i<<", "<<j<<", "<<k
                <<") is not stretchfree!\n"
                <<"|dR| - 1 = "<<R.determinant() - 1.<<endl
                <<"|R_OLD| - 1 = "<<RotationsOLD(i,j,k).determinant() - 1.<<endl
                <<"|R_NEW| - 1 = "<<OR.Quaternions(i,j,k).RotationMatrix.determinant() - 1.<<endl;

            cout<<"before correction:"<<endl<<"R:"<<endl<<R.print()<<endl;
            cout<<"correction factor: "<< pow(R.determinant(), 1./3.)<<endl;
            R /= pow(R.determinant(), 1./3.);
            cout<<"after correction"<<endl<<"R:"<<endl<<R.print()<<endl;
            cout<<"|dR|-1 = "<<R.determinant() - 1.<<endl;
//            exit(13);
        }
#endif
        PFP.PlasticStrain(i,j,k).rotate(R);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::ApplyHomExpansionToAccStress(
                                       ElasticProperties& EP, vStrain StrainInc)
{
    double invVolumeStrain = 1./(1+StrainInc[0] + StrainInc[1] + StrainInc[2]);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Stresses, 0,)
    {
        EP.Stresses(i,j,k) *= invVolumeStrain;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::SaveRotations(Orientations& OR,
                 BoundaryConditions& BC, Storage3D<dMatrix3x3, 0>& RotationsOLD)
{
    OR.SetBoundaryConditions(BC);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, OR.Quaternions, OR.Quaternions.Bcells(),)
    {
        RotationsOLD(i, j, k) = OR.Quaternions(i,j,k).RotationMatrix;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::SaveCRSS(PlasticFlowCP& PFCP, Storage3D<NodeVn<12>, 0>& CRSSOLD)
{
    int Nx = PFCP.Nx;
    int Ny = PFCP.Ny;
    int Nz = PFCP.Nz;

    if(CRSSOLD.IsNotAllocated())
        CRSSOLD.Allocate(Nx, Ny, Nz, 0);

    if(CRSSOLD.sizeX() != PFCP.CRSS.sizeX() or
        CRSSOLD.sizeY() != PFCP.CRSS.sizeY() or
        CRSSOLD.sizeZ() != PFCP.CRSS.sizeZ())
        CRSSOLD.Reallocate(Nx, Ny, Nz);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, PFCP.CRSS, 0,)
    {
        CRSSOLD(i,j,k) = PFCP.CRSS(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::SaveDamage(DamageProperties& DP,
                                                Storage3D<double, 0>& DamageOLD)
{
    int Nx = DP.Nx;
    int Ny = DP.Ny;
    int Nz = DP.Nz;

    if(DamageOLD.IsNotAllocated())
        DamageOLD.Allocate(Nx, Ny, Nz, 0);

    if(DamageOLD.sizeX() != DP.EffectiveDamage.sizeX() or
       DamageOLD.sizeY() != DP.EffectiveDamage.sizeY() or
       DamageOLD.sizeZ() != DP.EffectiveDamage.sizeZ())
        DamageOLD.Reallocate(Nx, Ny, Nz);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DP.EffectiveDamage, 0,)
    {
        DamageOLD(i,j,k) = DP.EffectiveDamage(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::SaveDamageN(DamageProperties& DP,
	Storage3D<double, 0>& DamageOLD, Storage3D<double, 0>& DamageNOLD)
{
	int Nx = DP.Nx;
	int Ny = DP.Ny;
	int Nz = DP.Nz;

	if (DamageOLD.IsNotAllocated())
		DamageOLD.Allocate(Nx, Ny, Nz, 0);

	if (DamageOLD.sizeX() != DP.EffectiveDamage.sizeX() or
		DamageOLD.sizeY() != DP.EffectiveDamage.sizeY() or
		DamageOLD.sizeZ() != DP.EffectiveDamage.sizeZ())
		DamageOLD.Reallocate(Nx, Ny, Nz);

	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DP.EffectiveDamage, 0, )
	{
		DamageOLD(i, j, k) = DP.EffectiveDamage(i, j, k);
	}
	OMP_PARALLEL_STORAGE_LOOP_END

	if (DamageNOLD.IsNotAllocated())
		DamageNOLD.Allocate(Nx, Ny, Nz, 0);

	if (DamageNOLD.sizeX() != DP.EffectiveDamage.sizeX() or
		DamageNOLD.sizeY() != DP.EffectiveDamage.sizeY() or
		DamageNOLD.sizeZ() != DP.EffectiveDamage.sizeZ())
		DamageNOLD.Reallocate(Nx, Ny, Nz);

	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DP.EffectiveDamageN, 0, )
	{
		DamageNOLD(i, j, k) = DP.EffectiveDamageN(i, j, k);
	}
	OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::SavePhaseField(PhaseField& Phase, PhaseField& PhaseOLD)
{
    PhaseOLD = Phase;
    PhaseOLD.thisclassname = "OLDPhaseField";
}

void LargeDeformationsMethods::SubtractPlasticStrainOLD(
              PlasticProperties& PFP, Storage3D<vStrain, 0>& PlasticStrainOLD)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, PFP.PlasticStrain, 0,)
    {
        PFP.PlasticStrain(i,j,k) -= PlasticStrainOLD(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::ADDPlasticStrainOLD(
              PlasticProperties& PFP, Storage3D<vStrain, 0>& PlasticStrainOLD)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, PFP.PlasticStrain, 0,)
    {
        PFP.PlasticStrain(i,j,k) += PlasticStrainOLD(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::DeletePlasticStrain(PlasticProperties& PFP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, PFP.PlasticStrain, 0,)
    {
        PFP.PlasticStrain(i,j,k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::RestoreCRSS(PlasticFlowCP& PFCP,
                                            Storage3D<NodeVn<12>, 0>& CRSSOLD)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, PFCP.CRSS, 0,)
    {
        PFCP.CRSS(i,j,k) = CRSSOLD(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::RestoreDamage(DamageProperties& DP,
                                            Storage3D<double, 0>& DamageOLD)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DP.EffectiveDamage, 0,)
    {
        DP.EffectiveDamage(i,j,k) = DamageOLD(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::RestoreDamageN(DamageProperties& DP,
	Storage3D<double, 0>& DamageOLD, Storage3D<double, 0>& DamageNOLD)
{
	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DP.EffectiveDamage, 0, )
	{
		DP.EffectiveDamage(i, j, k) = DamageOLD(i, j, k);
		DP.EffectiveDamageN(i, j, k) = DamageNOLD(i, j, k);
	}
	OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::RestoreRotations(Orientations& OR,
	BoundaryConditions& BC, Storage3D<dMatrix3x3, 0>& RotationsOLD)
{
	OR.SetBoundaryConditions(BC);
	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, OR.Quaternions, OR.Quaternions.Bcells(), )
	{
		OR.Quaternions(i, j, k).set(RotationsOLD(i, j, k));
	}
	OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::SaveEffectiveEigenStrains(PhaseField& Phase,
                                  ElasticProperties& EP, BoundaryConditions& BC,
                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD)
{
    EP.SetBoundaryConditions(BC);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields, Phase.Fields.Bcells(),)
    {
        EffectiveEigenStrainsOLD(i, j, k).set_to_zero();
        if(!(Phase.Interface(i, j, k)))
        {
            int alphaIdx = Phase.Fields(i,j,k).front().index;
            EffectiveEigenStrainsOLD(i, j, k) = EP.EigenStrains[alphaIdx];
        }
        else
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            {
                EffectiveEigenStrainsOLD(i, j, k) +=
                                     EP.EigenStrains[alpha->index]*alpha->value;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::SaveEffectiveEigenStrains(PhaseField& Phase,
                ElasticProperties& EP, Orientations& OR, BoundaryConditions& BC,
                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD)
{
    EP.SetBoundaryConditions(BC);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EffectiveEigenStrainsOLD, EffectiveEigenStrainsOLD.Bcells(),)
    {
        EffectiveEigenStrainsOLD(i, j, k).set_to_zero();
        if(!(Phase.Interface(i, j, k)))
        {
            int alphaIdx = Phase.Fields(i,j,k).front().index;
            EffectiveEigenStrainsOLD(i,j,k) = EP.EigenStrains[alphaIdx];
        }
        else
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            {
                EffectiveEigenStrainsOLD(i, j, k) +=
                                EP.EigenStrains[alpha->index]*alpha->value;
            }
        }
        EffectiveEigenStrainsOLD(i,j,k).rotate(OR.Quaternions(i,j,k).RotationMatrix);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::SaveEffectiveElasticConstants(
                  ElasticProperties& EP, Storage3D<dMatrix6x6, 0>& StiffnessOLD)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.EffectiveElasticConstants, EP.EffectiveElasticConstants.Bcells(),)
    {
        StiffnessOLD(i,j,k) = EP.EffectiveElasticConstants(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::SaveEffectiveEigenStrains(PhaseField& Phase,
                 Composition& Cx, ElasticProperties& EP, BoundaryConditions& BC,
                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD)
{
    int Ncomp = EP.Ncomp;

    EP.SetBoundaryConditions(BC);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields, Phase.Fields.Bcells(),)
    {
        EffectiveEigenStrainsOLD(i,j,k).set_to_zero();
        if(!(Phase.Interface(i,j,k)))
        {
            int index = Phase.Fields(i,j,k).front().index;;
            int pIndex = Phase.FieldsStatistics[index].Phase;

            EffectiveEigenStrainsOLD(i,j,k) = EP.EigenStrains[index];

            for(int comp = 0; comp < Ncomp; comp++)
            {
                double delta = (Cx.Phase(i,j,k)({pIndex, comp}) -
                                                   EP.Cref({pIndex, comp}));
                for(int dir = 0; dir < 6; dir++)
                {
                    EffectiveEigenStrainsOLD(i,j,k)[dir] += delta *
                                              EP.Lambda({index, comp})[dir];
                }
            }
        }
        else
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            {
                int index = alpha->index;
                int pIndex = Phase.FieldsStatistics[index].Phase;
                double phiAlpha = alpha->value;

                EffectiveEigenStrainsOLD(i, j, k) +=
                                                EP.EigenStrains[index]*phiAlpha;

                for(int comp = 0; comp < Ncomp; comp++)
                {
                    double delta = (Cx.Phase(i,j,k)({pIndex, comp}) -
                                                         EP.Cref({pIndex, comp}));

                    for(int dir = 0; dir < 6; dir++)
                    {
                        EffectiveEigenStrainsOLD(i, j, k)[dir] += delta*
                                       phiAlpha*EP.Lambda({index, comp})[dir];
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::SaveEffectiveKirkendallEigenStrains(PhaseField& Phase,
                 Composition& Cx, ElasticProperties& EP, BoundaryConditions& BC,
                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD)
{
    int Ncomp = EP.Ncomp;

    EP.SetBoundaryConditions(BC);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields, Phase.Fields.Bcells(),)
    {
        EffectiveEigenStrainsOLD(i, j, k).set_to_zero();
        if(!(Phase.Interface(i, j, k)))
        {
            int index = Phase.Fields(i,j,k).front().index;;
            int pIndex = Phase.FieldsStatistics[index].Phase;

            EffectiveEigenStrainsOLD(i, j, k) = EP.EigenStrains[index];

            double VolumeStrain = -1.0;
            for(int comp = 0; comp < Ncomp; comp++)
            {
                VolumeStrain += Cx.Phase(i,j,k)({pIndex, comp})*
                        Cx.MolarVolume({pIndex, comp});
            }

            double KirkendallStrain = VolumeStrain*0.3333333333333333333333;
            for(int dir = 0; dir < 3; dir++)
            {
                EffectiveEigenStrainsOLD(i,j,k)[dir] += KirkendallStrain;
            }
        }
        else
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            {
                int index = alpha->index;
                int pIndex = Phase.FieldsStatistics[index].Phase;
                double phiAlpha = alpha->value;

                EffectiveEigenStrainsOLD(i, j, k) +=
                                                EP.EigenStrains[index]*phiAlpha;

                double VolumeStrain = -1.0;
                for(int comp = 0; comp < Ncomp; comp++)
                {
                    VolumeStrain += Cx.Phase(i,j,k)({pIndex, comp})*
                                    Cx.MolarVolume({pIndex, comp});
                }

                double KirkendallStrain = VolumeStrain*0.3333333333333333333333;
                for(int dir = 0; dir < 3; dir++)
                {
                    EffectiveEigenStrainsOLD(i,j,k)[dir] +=
                                                  KirkendallStrain*phiAlpha;
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

double LargeDeformationsMethods::GetEigenstrainNorm(std::vector<OPObject*> Objects,
                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD,
                                         Storage3D<dMatrix6x6, 0>& StiffnessOLD,
                                    double dt, int nSubsteps, double MaxAllowedStrainIncrement,
                                    bool verbose, bool screenWrite)
{
    PhaseField* Phase           = (static_cast<PhaseField*>(OPObject::findOPObject(Objects, "PhaseField", "LargeDeformationsMethods", "GetEigenstrainNorm", true, verbose)));
    ElasticProperties* EP       = (static_cast<ElasticProperties*>(OPObject::findOPObject(Objects, "ElasticProperties", "LargeDeformationsMethods", "GetEigenstrainNorm", true, verbose)));
    Orientations* OR            = (static_cast<Orientations*>(OPObject::findOPObject(Objects, "Orientations", "LargeDeformationsMethods", "Solve", true, verbose)));

    PlasticProperties* PFP    = (static_cast<PlasticProperties*>(OPObject::findOPObject(Objects, "PlasticProperties", "LargeDeformationsMethods", "GetEigenstrainNorm", false, verbose)));

    double StrainNorm = 1.0;

    double GlobalMaxStrainIncrement = 0.0;

    #pragma omp parallel //OMP BEGIN
    {
        double MaxStrainIncrement = 0.0;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP->EffectiveEigenStrains, 0,)
        {
            EP->EffectiveEigenStrains(i, j, k).set_to_zero();
            if(!(Phase->Interface(i, j, k)))
            {
				int alphaIdx = Phase->Fields(i, j, k).front().index;
				vStrain temp = EP->EigenStrains[alphaIdx].rotated(OR->Quaternions(i, j, k).RotationMatrix).Ln();
				vStrain temp2 = EffectiveEigenStrainsOLD(i, j, k).Ln();
                
				EP->EffectiveEigenStrains(i, j, k) = temp - temp2;
                /*EP->EffectiveEigenStrains(i, j, k) =
                            EP->EigenStrains[alphaIdx].rotated(OR->Quaternions(i,j,k).RotationMatrix).Ln() -
                            EffectiveEigenStrainsOLD(i,j,k).Ln();*/
            }
            else
            {
                for(auto alpha = Phase->Fields(i, j, k).cbegin();
                         alpha != Phase->Fields(i,j,k).cend(); ++alpha)
                {
                    EP->EffectiveEigenStrains(i, j, k) +=
                             EP->EigenStrains[alpha->index].rotated(OR->Quaternions(i,j,k).RotationMatrix)*alpha->value;
                }
                EP->EffectiveEigenStrains(i,j,k) =
                                    EP->EffectiveEigenStrains(i,j,k).Ln() -
                                    EffectiveEigenStrainsOLD(i,j,k).Ln();
            }
            if(PFP != nullptr)
                EP->EffectiveEigenStrains(i,j,k) += PFP->PlasticStrain(i,j,k);    /// plastic strain increment is already logarithmic

            if(StiffnessOLD(i,j,k)(0,0)>0.0)
            {
                dMatrix6x6 dS = EP->EffectiveElasticConstants(i,j,k).inverted()-
                                StiffnessOLD(i,j,k).inverted();
                EP->EffectiveEigenStrains(i,j,k) +=
                                            (dS*EP->Stresses(i,j,k)).Ln();
            }

            for(int dir = 0; dir < 6; dir++)
            {
                double StrainIncrement = fabs(EP->EffectiveEigenStrains(i, j, k)[dir] / nSubsteps);
                if(StrainIncrement > MaxStrainIncrement)
                    MaxStrainIncrement = StrainIncrement;
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        #pragma omp critical
        if(MaxStrainIncrement > GlobalMaxStrainIncrement)
            GlobalMaxStrainIncrement = MaxStrainIncrement;
    }

    if(GlobalMaxStrainIncrement > MaxAllowedStrainIncrement)
    {
        int iterations = ceil(GlobalMaxStrainIncrement/
                              MaxAllowedStrainIncrement);
        StrainNorm = 1./iterations;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP->EffectiveEigenStrains, 0,)
        {
                EP->EffectiveEigenStrains(i,j,k) *= StrainNorm;
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    stringstream message;
    message << GlobalMaxStrainIncrement << " (" << std::setprecision(3)
            << GlobalMaxStrainIncrement/MaxAllowedStrainIncrement << "fold of max. allowed inc.)";

    if (screenWrite) Info::WriteStandard("GlobalMaxStrainIncrement", message.str());
    return StrainNorm;
}

double LargeDeformationsMethods::SetEigenstrainsIncrement(std::vector<OPObject*> Objects,
                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD,
                                         Storage3D<dMatrix6x6, 0>& StiffnessOLD,
                                                        double dt, int nSubsteps,
                                                        bool verbose)
{
    PhaseField* Phase           = (static_cast<PhaseField*>(OPObject::findOPObject(Objects, "PhaseField", "LargeDeformationsMethods", "SetEigenstrainsIncrement", true, verbose)));
    ElasticProperties* EP       = (static_cast<ElasticProperties*>(OPObject::findOPObject(Objects, "ElasticProperties", "LargeDeformationsMethods", "SetEigenstrainsIncrement", true, verbose)));
    Orientations* OR            = (static_cast<Orientations*>(OPObject::findOPObject(Objects, "Orientations", "LargeDeformationsMethods", "SetEigenstrainsIncrement", true, verbose)));

    PlasticProperties* PFP    = (static_cast<PlasticProperties*>(OPObject::findOPObject(Objects, "PlasticProperties", "LargeDeformationsMethods", "SetEigenstrainNorm", false, verbose)));

    BoundaryConditions* BC      = (static_cast<BoundaryConditions*>(OPObject::findOPObject(Objects, "BoundaryConditions", "SetEigenstrainsIncrement", "SetEigenstrainNorm", true, verbose)));



    double StrainNorm = 1.0;

    StrainNorm /= nSubsteps;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP->EffectiveEigenStrains, 0,)
    {
        EP->EffectiveEigenStrains(i, j, k).set_to_zero();
        if(!(Phase->Interface(i, j, k)))
        {
            int alphaIdx = Phase->Fields(i,j,k).front().index;
            EP->EffectiveEigenStrains(i, j, k) =
                                (EP->EigenStrains[alphaIdx].rotated(OR->Quaternions(i,j,k).RotationMatrix).Ln() -
                                EffectiveEigenStrainsOLD(i,j,k).Ln()) * StrainNorm;
        }
        else
        {
            for(auto alpha = Phase->Fields(i, j, k).cbegin();
                     alpha != Phase->Fields(i,j,k).cend(); ++alpha)
            {
                EP->EffectiveEigenStrains(i, j, k) +=
                             EP->EigenStrains[alpha->index].rotated(OR->Quaternions(i,j,k).RotationMatrix) * alpha->value;
            }
            EP->EffectiveEigenStrains(i,j,k) =
                                (EP->EffectiveEigenStrains(i,j,k).Ln() -
                                EffectiveEigenStrainsOLD(i,j,k).Ln()) * StrainNorm;
        }
        if(StiffnessOLD(i,j,k)(0,0) > 0.0)
        {
            dMatrix6x6 dS = (EP->EffectiveElasticConstants(i,j,k).inverted() -
                            StiffnessOLD(i,j,k).inverted()) * StrainNorm;
            EP->EffectiveEigenStrains(i,j,k) +=
                                    (dS*EP->Stresses(i,j,k)).Ln();
        }
        //EP->EffectiveEigenStrains(i, j, k) *= StrainNorm;
		if (PFP != nullptr)
		{
			EP->EffectiveEigenStrains(i, j, k) += PFP->PlasticStrain(i, j, k);        /// plastic strain increment is already logarithmic and correspond to the reduced timestep
		}
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    //EP->SmoothEigenstrain();
    EP->SetBoundaryConditions(*BC);
    return StrainNorm;
}

double LargeDeformationsMethods::GetEigenstrainNorm(PhaseField& Phase,
                                         Composition& Cx, ElasticProperties& EP,
                                Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD,
                                               double MaxAllowedStrainIncrement)
{
    int Ncomp = EP.Ncomp;

    double GlobalMaxStrainIncrement = 0.0;

    double StrainNorm = 1.0;

    #pragma omp parallel //OMP BEGIN
    {
        double MaxStrainIncrement = 0.0;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields, Phase.Fields.Bcells(),)
        {
            vStrain LocEigenStrain;
            LocEigenStrain.set_to_zero();

            if(!(Phase.Interface(i, j, k)))
            {
                int index = Phase.Fields(i,j,k).front().index;;
                int pIndex = Phase.FieldsStatistics[index].Phase;
                LocEigenStrain = EP.EigenStrains[index];

                for(int comp = 0; comp < Ncomp; comp++)
                {
                    double delta = (Cx.Phase(i,j,k)({pIndex, comp}) -
                                                       EP.Cref({pIndex, comp}));
                    for(int dir = 0; dir < 6; dir++)
                    {
                        LocEigenStrain[dir] += delta*
                                                    EP.Lambda({index, comp})[dir];
                    }
                }
            }
            else
            {
                for(auto alpha = Phase.Fields(i, j, k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend(); ++alpha)
                {
                    int index = alpha->index;
                    int pIndex = Phase.FieldsStatistics[index].Phase;
                    double phiAlpha = alpha->value;

                    LocEigenStrain += EP.EigenStrains[index]*phiAlpha;

                    for(int comp = 0; comp < Ncomp; comp++)
                    {
                        double delta = (Cx.Phase(i,j,k)({pIndex, comp}) -
                                                       EP.Cref({pIndex, comp}));

                        for(int dir = 0; dir < 6; dir++)
                        {
                            LocEigenStrain[dir] += delta*phiAlpha*
                                                    EP.Lambda({index, comp})[dir];
                        }
                    }
                }
            }
            LocEigenStrain = LocEigenStrain.Ln() -
                                           EffectiveEigenStrainsOLD(i,j,k).Ln();

        }
        OMP_PARALLEL_STORAGE_LOOP_END

        #pragma omp critical
        if(MaxStrainIncrement > GlobalMaxStrainIncrement)
        {
            GlobalMaxStrainIncrement = MaxStrainIncrement;
        }
    }//OMP END

    if(GlobalMaxStrainIncrement > MaxAllowedStrainIncrement)
    {
        int iterations = ceil(GlobalMaxStrainIncrement/MaxAllowedStrainIncrement);
        StrainNorm = 1.0/iterations;
    }

    return StrainNorm;
}

double LargeDeformationsMethods::SetEigenstrainsIncrement(PhaseField& Phase,
                                         Composition& Cx, ElasticProperties& EP,
                  Storage3D<vStrain, 0>& EffectiveEigenStrainsOLD, int substeps)
{
    int Ncomp = EP.Ncomp;

    double StrainNorm = 1.0;

    StrainNorm /= substeps;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields, Phase.Fields.Bcells(),)
    {
        EP.EffectiveEigenStrains(i, j, k).set_to_zero();
        if(!(Phase.Interface(i, j, k)))
        {
            int index = Phase.Fields(i,j,k).front().index;;
            int pIndex = Phase.FieldsStatistics[index].Phase;
            EP.EffectiveEigenStrains(i, j, k) = EP.EigenStrains[index];

            for(int comp = 0; comp < Ncomp; comp++)
            {
                double delta = (Cx.Phase(i,j,k)({pIndex, comp}) -
                                                 EP.Cref({pIndex, comp}));
                for(int dir = 0; dir < 6; dir++)
                {
                    EP.EffectiveEigenStrains(i, j, k)[dir] += delta *
                                            EP.Lambda({index, comp})[dir];
                }
            }
        }
        else
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend(); ++alpha)
            {
                int index = alpha->index;
                int pIndex = Phase.FieldsStatistics[index].Phase;
                double phiAlpha = alpha->value;

                EP.EffectiveEigenStrains(i, j, k) +=
                                        EP.EigenStrains[index]*phiAlpha;

                for(int comp = 0; comp < Ncomp; comp++)
                {
                    double delta = (Cx.Phase(i,j,k)({pIndex, comp}) -
                                                 EP.Cref({pIndex, comp}));

                    for(int dir = 0; dir < 6; dir++)
                    {
                        EP.EffectiveEigenStrains(i, j, k)[dir] += delta*
                                   phiAlpha*EP.Lambda({index, comp})[dir];
                    }
                }
            }
        }
        EP.EffectiveEigenStrains(i,j,k) =
                   (EP.EffectiveEigenStrains(i,j,k).Ln() -
                    EffectiveEigenStrainsOLD(i,j,k).Ln())*StrainNorm;

        EP.EffectiveEigenStrains(i, j, k) *= StrainNorm;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    return StrainNorm;
}

double LargeDeformationsMethods::GetAppliedStrainNorm(ElasticProperties& EP,
                                               double MaxAllowedStrainIncrement,
                                               bool screenWrite)
{
    /*
     * Calculates the proper number of substeps and returns the norm factor, the
     * EP.AppliedStrain and EP.Eigenstrains are to be multiplied with.
     */
    double StrainNorm = 1.0;
    double MaxStrainIncrement = 0.0;
	vStrain StrainrateInc = EP.AppliedStrainRate * EP.dt;
    EP.EffectiveAppliedStrain = (EP.AppliedStrain.Ln() + StrainrateInc.Ln()) -
                                EP.AppliedStrainOLD.Ln();

    for(int dir = 0; dir < 6; dir++)
    {
        if(fabs(EP.EffectiveAppliedStrain[dir]) > MaxStrainIncrement)
                  MaxStrainIncrement = fabs(EP.EffectiveAppliedStrain[dir]);
    }

    if(MaxStrainIncrement > MaxAllowedStrainIncrement)
    {
        int iterations = ceil(MaxStrainIncrement/MaxAllowedStrainIncrement);
        StrainNorm = 1.0/iterations;

        EP.EffectiveAppliedStrain *= StrainNorm;
    }

    stringstream message;
    message << MaxStrainIncrement << " (" << std::setprecision(3)
            << MaxStrainIncrement/MaxAllowedStrainIncrement << "fold of max. allowed inc.)";
    if (screenWrite) Info::WriteStandard("MaxAppliedStrainIncrement", message.str());

    return StrainNorm;
}

double LargeDeformationsMethods::GetAppliedStressNorm(ElasticProperties& EP,
                                               double MaxAllowedStrainIncrement)
{
    /*
     * Calculates the proper number of substeps and returns the norm factor, the
     * EP.AppliedStrain and EP.Eigenstrains are to be multiplied with.
     */

    double minC = 1.0e30; double maxC = 1.0;
    dMatrix6x6 minCx;
    dMatrix6x6 maxCx;

    for(int pf = 0; pf < EP.ElasticConstants.Size(); ++pf)
    {
        if (EP.ElasticConstants[pf].norm() > maxC)
        {
            maxC = EP.ElasticConstants[pf].norm();
            maxCx = EP.ElasticConstants[pf];
        }
        if (EP.ElasticConstants[pf].norm() < minC and EP.ElasticConstants[pf].norm() > 0)
        {
            minC = EP.ElasticConstants[pf].norm();
            minCx = EP.ElasticConstants[pf];
        }
    }
    dMatrix6x6 C0 = (maxCx + minCx)*0.5;

    vStrain approxAppliedStrain = C0.inverted()*(EP.AppliedStress - EP.AppliedStressOLD);

//    cout << "AppStress: " << EP.AppliedStress.print() << endl;
//    cout << "AppStressOLD: " << EP.AppliedStressOLD.print() << endl;
//    cout << "Guessed: " << approxAppliedStrain.print() << endl;

    double StressNorm = 1.0;
    double MaxStrainIncrement = 0.0;
    EP.EffectiveAppliedStrain = approxAppliedStrain/*.Ln() -
                                EP.AppliedStrainOLD.Ln()*/;

    for(int dir = 0; dir < 6; dir++)
    {
        if(fabs(EP.EffectiveAppliedStrain[dir]) > MaxStrainIncrement)
                  MaxStrainIncrement = fabs(EP.EffectiveAppliedStrain[dir]);
    }

    EP.EffectiveAppliedStrain.set_to_zero();

    if(MaxStrainIncrement > MaxAllowedStrainIncrement)
    {
        int iterations = ceil(MaxStrainIncrement/MaxAllowedStrainIncrement);
        StressNorm = 1.0/iterations;

        EP.EffectiveAppliedStress = EP.AppliedStress*StressNorm;
    }

    stringstream message;
    message << MaxStrainIncrement << " (" << std::setprecision(3)
            << MaxStrainIncrement/MaxAllowedStrainIncrement << "fold of max. allowed inc.)";
    Info::WriteStandard("MaxAppliedStrainIncrement (StressBC)", message.str());

    return StressNorm;
}

double LargeDeformationsMethods::SetAppliedStrainIncrement(ElasticProperties& EP,
                                                                   int substeps)
{
    /*
     * Sets the EP.EffectiveAppliedStrain to EP.AppliedStrain/substeps.
     */
    double StrainNorm = 1.0;
    StrainNorm /= substeps;
	vStrain temp = EP.AppliedStrain + EP.AppliedStrainRate*EP.dt;
    EP.EffectiveAppliedStrain = (temp.Ln() -
                                 EP.AppliedStrainOLD.Ln())*StrainNorm;
	EP.AppliedStrain = temp;
	//cout << "EffectiveAppliedStrain = " << EP.EffectiveAppliedStrain.print()  << " = ( " << EP.AppliedStrain.Ln().print() << " - "
	//	<< EP.AppliedStrainOLD.Ln().print()<<" )* " << StrainNorm<< endl;

    return StrainNorm;
}

double LargeDeformationsMethods::SetAppliedStressIncrement(ElasticProperties& EP,
                                                                   int substeps)
{
    /*
     * Sets the EP.EffectiveAppliedStrain to EP.AppliedStrain/substeps.
     */
    double StrainNorm = 1.0;
    StrainNorm /= substeps;
    EP.EffectiveAppliedStress = (EP.AppliedStress - EP.AppliedStressOLD)*StrainNorm;

    return StrainNorm;
}

void LargeDeformationsMethods::SetHomogeneousVelocity(PhaseField& Phase,
                          Velocities& Vel, vStrain HomogeneousStrain, double dt)
{
    int Nx = Vel.Nx;
    int Ny = Vel.Ny;
    int Nz = Vel.Nz;

    double dx = Vel.dx;
    double dtInv = 1.0/dt;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, Vel.Average, Vel.Average.Bcells(),)
    {
        Vel.Average(i,j,k)[0] = HomogeneousStrain[0]*(i - (Nx+1.0)*0.5)*dx*dtInv;
        Vel.Average(i,j,k)[1] = HomogeneousStrain[1]*(j - (Ny+1.0)*0.5)*dx*dtInv;
        Vel.Average(i,j,k)[2] = HomogeneousStrain[2]*(k - (Nz+1.0)*0.5)*dx*dtInv;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::SetHomogeneousVelocityBC(Velocities& Vel)
{
    int Nx = Vel.Nx;
    int Ny = Vel.Ny;
    int Nz = Vel.Nz;

    int Bcells = Vel.Average.Bcells();

    for(int i = -Bcells; i < Nx+Bcells; i++)
    for(int k = -Bcells; k < Nz+Bcells; k++)
    {
        Vel.Average(i,     -Bcells, k) = Vel.Average(i, 0, k)*2.0 - Vel.Average(i, 1, k);
        Vel.Average(i, Ny+Bcells-1, k) = Vel.Average(i, Ny, k)*2.0 - Vel.Average(i, Ny-1, k);
    }

    for(int j = -Bcells; j < Ny+Bcells; j++)
    for(int k = -Bcells; k < Nz+Bcells; k++)
    {
        Vel.Average(    -Bcells, j, k) = Vel.Average( 0, j, k)*2.0 - Vel.Average(   1, j, k);
        Vel.Average(Nx+Bcells-1, j, k) = Vel.Average(Nx, j, k)*2.0 - Vel.Average(Nx-1, j, k);
    }

    for(int i = -Bcells; i < Nx+Bcells; i++)
    for(int j = -Bcells; j < Ny+Bcells; j++)
    {
        Vel.Average(i, j,     -Bcells) = Vel.Average(i, j,  0)*2.0 - Vel.Average(i, j,    1);
        Vel.Average(i, j, Nz+Bcells-1) = Vel.Average(i, j, Nz)*2.0 - Vel.Average(i, j, Nz-1);
    }
}

void LargeDeformationsMethods::SetEffectiveElasticConstants(PhaseField& Phase, ElasticProperties& EP, Orientations& OR)
{
    EP.SetEffectiveElasticConstants(Phase);
    EP.ApplyRotation(OR);

    const int TensorIdx[6][2] = {{0, 0}, {1, 1}, {2, 2}, {1, 2}, {0, 2}, {0, 1}};

    const int Nx = EP.Nx;
    const int Ny = EP.Ny;
    const int Nz = EP.Nz;

    for(int i = 0; i < Nx; i++)
    for(int j = 0; j < Ny; j++)
    for(int k = 0; k < Nz; k++)
    {
        const dMatrix3x3 sigma = EP.Stresses(i,j,k).tensor();

        for(int m = 0; m < 6; m++)
        {
            for(int n = 0; n < 6; n++)
            {
                const int ii = TensorIdx[m][0];
                const int jj = TensorIdx[m][1];
                const int kk = TensorIdx[n][0];
                const int ll = TensorIdx[n][1];
                EP.EffectiveElasticConstants(i,j,k)(m,n) +=
                                      sigma(ii,kk)*(jj==ll) + sigma(jj,kk)*(ii==ll)
                                    + sigma(ii,ll)*(jj==kk) + sigma(jj,ll)*(ii==kk);
                if(EP.EffectiveElasticConstants(i,j,k)(m,n) < -1.0e8)
                                EP.EffectiveElasticConstants(i,j,k)(m,m) = -1.0e8;
            }
            if(EP.EffectiveElasticConstants(i,j,k)(m,m) < 1.0e6)
               EP.EffectiveElasticConstants(i,j,k)(m,m) = 1.0e6;
        }
    }
}

void LargeDeformationsMethods::RestoreInterface(PhaseField& Phase,
                                 BoundaryConditions& BC, InterfaceEnergy& Sigma,
                          InterfaceMobility& Mu, InterfaceField& Psi, double dt)
{
    Sigma.CalculateCubic(Phase);
    Mu.CalculateCubic(Phase);
    int Nx = Phase.Nx;
    int Ny = Phase.Ny;
    int Nz = Phase.Nz;
    double Prefactor2 = Pi*Pi/(Phase.Eta*Phase.Eta);
    #pragma omp parallel //OMP BEGIN
    {
        int nThreads = omp_get_num_threads();
        int myThread = omp_get_thread_num();

        int xBeg = (myThread*Nx)/nThreads;
        int xEnd = ((myThread+1)*Nx)/nThreads;

        for(int i = xBeg; i < xEnd; i++)
        for(int j =    0; j < Ny; j++)
        for(int k =    0; k < Nz; k++)
        if (Phase.Fields(i,j,k).flag)
        {
            double sum = 1.0;
            for (auto alpha = Phase.Fields(i,j,k).cbegin();
                      alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
            for (auto  beta = alpha + 1;
                       beta < Phase.Fields(i,j,k).cend(); ++beta)
            {
                sum -= alpha->value*beta->value;
            }

            if(sum > 0.0 and sum< 1.0)
            {
                double norm_1 = 1.0/double(Phase.Fields(i,j,k).size());
                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta != Phase.Fields(i,j,k).cend(); ++beta)
                {
                    double dPsi_dt  = Sigma(i,j,k, alpha->index, beta->index)*
                                      ((alpha->laplacian + Prefactor2*alpha->value) -
                                       ( beta->laplacian + Prefactor2* beta->value));
                    if(Phase.Fields(i,j,k).size() > 2)
                    for(auto gamma = Phase.Fields(i,j,k).cbegin();
                             gamma != Phase.Fields(i,j,k).cend(); ++gamma)
                    if((gamma != alpha) && (gamma != beta))
                    {
                        dPsi_dt += (Sigma(i,j,k,  beta->index, gamma->index) -
                                    Sigma(i,j,k, alpha->index, gamma->index))*
                                    (gamma->laplacian + Prefactor2*gamma->value);
                    }
                    dPsi_dt *= Mu(i,j,k, alpha->index, beta->index)*norm_1;
                    Phase.FieldsDot(i,j,k).add_asym(alpha->index, beta->index,  dPsi_dt);
                }
            }
        }
    }//OMP END
    Phase.NormalizeIncrements(BC, dt);
    Phase.MergeIncrements(BC, dt);
}

void LargeDeformationsMethods::SetPlasticStrainEquivalent(
                                  Storage3D<double, 0>& PlasticStrainEquivalent,
                                           Storage3D<vStrain, 0>& PlasticStrain)
{
    const double factor = sqrt(2.0/3.0);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, PlasticStrain, 0,)
    {
        PlasticStrainEquivalent(i,j,k) = factor*PlasticStrain(i,j,k).norm();
        if(factor*PlasticStrain(i,j,k).norm() < 0.0)
            cout<<"fehler!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::SetPlasticStrainEquivalentIntf(PhaseField& Phase,
                                  Storage3D<double, 0>& PlasticStrainEquivalent,
                           Storage3D<vStrain, 0>& PlasticStrain, int phaseIndex)
{
    const double factor = sqrt(2.0/3.0);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, PlasticStrain, 0,)
    if(Phase.Fields(i,j,k).get(phaseIndex) > 0.0)
    {
        PlasticStrainEquivalent(i,j,k) = factor*PlasticStrain(i,j,k).norm();
        if(factor*PlasticStrain(i,j,k).norm() < 0.0)
            cout<<"fehler!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::CalculatePhaseDamageEnergy(
                   std::vector<double>& PhaseDamageEnergy, DamageProperties& DP,
                                       ElasticProperties& EP, PhaseField& Phase,
                                                Storage3D<double, 0>& DamageOLD)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.EffectiveElasticConstants, 0,)
    if(DP.EffectiveDamage(i,j,k) > 0.0001 and
       EP.EffectiveElasticConstants(i,j,k)(0,0) > 0.0 and
       DP.EffectiveDamage(i,j,k) - DamageOLD(i,j,k) > 0.0)
    {
        dMatrix6x6 StiffnesDecr = EP.EffectiveElasticConstants(i,j,k)*
                                   (DP.EffectiveDamage(i,j,k)-DamageOLD(i,j,k));
        vStrain locStrain = StiffnesDecr.inverted()*
                                                  EP.Stresses(i,j,k);
        double energy = locStrain*EP.Stresses(i,j,k);
        for (auto it = Phase.Fields(i,j,k).cbegin();
                  it != Phase.Fields(i,j,k).cend(); ++it)
        {
            PhaseDamageEnergy[it->index] += 0.5*it->value*energy/
                                       Phase.FieldsStatistics[it->index].Volume;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformationsMethods::WriteEffectiveElasticConstantsDiffVTK(ElasticProperties& EP, int tStep, Storage3D<dMatrix6x6, 0>& StiffnessOld)
{
    int Nx = EP.Nx;
    int Ny = EP.Ny;
    int Nz = EP.Nz;

    stringstream buffer;
    std::vector<int> DataTypes {1};
//    vector<long int> dimV {EP.EffectiveEigenStrains.sizeX(),
//        EP.EffectiveEigenStrains.sizeY(), EP.EffectiveEigenStrains.sizeZ()};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    for (auto m = 0; m < 6; m++)
    for (auto n = m; n < 6; n++)
    {
        buffer << "<DataArray type = \"Float64\" Name = \""
               << "C_" << m + 1 << n + 1
               << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < Nz; k++)
        for(int j = 0; j < Ny; j++)
        for(int i = 0; i < Nx; i++)
        {
            buffer << (EP.EffectiveElasticConstants(i,j,k)(m,n) - StiffnessOld(i,j,k)(m,n))<< endl;
        }
        buffer << "</DataArray>" << endl;
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "EffectiveElasticConstantsDiff", tStep);
}
} // namespace opensim

