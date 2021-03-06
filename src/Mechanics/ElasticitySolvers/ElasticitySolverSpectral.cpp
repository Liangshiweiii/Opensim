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

#include "fftw3.h"
#include "Settings.h"
#include "Info.h"
#include "Mechanics/ElasticitySolvers/ElasticitySolverSpectral.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "Orientations.h"
#include "Velocities.h"
#include "BoundaryConditions.h"

namespace opensim
{
using namespace std;

ElasticitySolverSpectral::ElasticitySolverSpectral(Settings& locSettings)
{
    this->Initialize(locSettings);
}

void ElasticitySolverSpectral::Initialize(Settings& locSettings)
{
    thisclassname = "ElasticitySolverSpectral";

    Nx  =  locSettings.Nx;
    Ny  =  locSettings.Ny;
    Nz  =  locSettings.Nz;
    Nz2 = (locSettings.Nz)/2+1;
    rlSIZE = Nx*Ny*Nz;
    rcSIZE = Nx*Ny*Nz2;

    dx = locSettings.dx;

    DPi_Nx = 2.0*Pi/double(Nx);
    DPi_Ny = 2.0*Pi/double(Ny);
    DPi_Nz = 2.0*Pi/double(Nz);

    Norm = 1.0/double(rlSIZE);
    // Arrays allocation:
    for(int n = 0; n < 6; n++)
    {
        rlRHSide[n]  = new double[rlSIZE] ();
        rcRHSide[n]  = new complex<double> [rcSIZE] ();
    }
    for(int n = 0; n < 3; n++)
    {
    	rlU[n] = new double[rlSIZE] ();
    	rcU[n] = new complex<double> [rcSIZE] ();
        Q[n]   = new double[rcSIZE] ();
    }
    for(int n = 0; n < 9; n++)
    {
        rlDefGrad[n]  = new double[rlSIZE] ();
        rcDefGrad[n]  = new complex<double> [rcSIZE] ();
    }

    SetQ();

//    removed since it slows down the initialisation!
//    fftw_init_threads();
//    fftw_plan_with_nthreads(omp_get_max_threads());

    for(int n = 0; n < 6; n++)
    {
        ForwardPlanRHS[n]  = fftw_plan_dft_r2c_3d
                        (Nx, Ny, Nz, rlRHSide[n],
                         reinterpret_cast<fftw_complex*> (rcRHSide[n]),
                         FFTW_PATIENT);
    }

    for(int n = 0; n < 9; n++)
    {
        BackwardPlanDefGrad[n] = fftw_plan_dft_c2r_3d
                        (Nx, Ny, Nz,
                         reinterpret_cast<fftw_complex*> (rcDefGrad[n]),
                         rlDefGrad[n],
                         FFTW_PATIENT);
    }

    for(int n = 0; n < 3; n++)
	{
		BackwardPlanU[n] = fftw_plan_dft_c2r_3d
						(Nx, Ny, Nz,
						 reinterpret_cast<fftw_complex*> (rcU[n]),
						 rlU[n],
						 FFTW_PATIENT);
	}

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void ElasticitySolverSpectral::Initialize(Settings& locSettings, BoundaryConditions& BC)
{
    thisclassname = "ElasticitySolverSpectral";

    if(BC.BC0X != Periodic or BC.BCNX != Periodic)
    {
        Nx  =  locSettings.Nx*2;
    }
    else
    {
        Nx  =  locSettings.Nx;
    }
    if(BC.BC0Y != Periodic or BC.BCNY != Periodic)
    {
        Ny  =  locSettings.Ny*2;
    }
    else
    {
        Ny  =  locSettings.Ny;
    }
    if(BC.BC0Z != Periodic or BC.BCNZ != Periodic)
    {
        Nz  =  locSettings.Nz*2;
    }
    else
    {
        Nz  =  locSettings.Nz;
    }

    Nz2 = Nz/2+1;
    rlSIZE = Nx*Ny*Nz;
    rcSIZE = Nx*Ny*Nz2;

    dx = locSettings.dx;

    DPi_Nx = 2.0*Pi/double(Nx);
    DPi_Ny = 2.0*Pi/double(Ny);
    DPi_Nz = 2.0*Pi/double(Nz);

    Norm = 1.0/double(rlSIZE);
    // Arrays allocation:
    for(int n = 0; n < 6; n++)
	{
		rlRHSide[n]  = new double[rlSIZE] ();
		memset(rlRHSide[n], 0, rlSIZE*sizeof(double));

		rcRHSide[n]  = new complex<double> [rcSIZE] ();
		memset(rcRHSide[n], 0, rcSIZE*sizeof(complex<double>));
	}
	for(int n = 0; n < 3; n++)
	{
		rlU[n] = new double[rlSIZE] ();
		memset(rlU[n], 0, rlSIZE*sizeof(double));

		rcU[n] = new complex<double> [rcSIZE] ();
		memset(rcU[n], 0, rcSIZE*sizeof(complex<double>));

		Q[n]   = new double[rcSIZE] ();
		memset(Q[n], 0, rcSIZE*sizeof(double));
	}
	for(int n = 0; n < 9; n++)
	{
		rlDefGrad[n] = new double[rlSIZE] ();
		memset(rlDefGrad[n], 0, rlSIZE*sizeof(double));

		rcDefGrad[n] = new complex<double> [rcSIZE] ();
		memset(rcDefGrad[n], 0, rcSIZE*sizeof(complex<double>));
	}

    SetQ();

//    removed since it slows down the initialisation!
//    fftw_init_threads();
//    fftw_plan_with_nthreads(omp_get_max_threads());

    for(int n = 0; n < 6; n++)
    {
        ForwardPlanRHS[n]  = fftw_plan_dft_r2c_3d
                        (Nx, Ny, Nz, rlRHSide[n],
                         reinterpret_cast<fftw_complex*> (rcRHSide[n]),
                         FFTW_PATIENT);
    }

    for(int n = 0; n < 9; n++)
    {
        BackwardPlanDefGrad[n] = fftw_plan_dft_c2r_3d
                        (Nx, Ny, Nz,
                         reinterpret_cast<fftw_complex*> (rcDefGrad[n]),
                         rlDefGrad[n],
                         FFTW_PATIENT);
    }

    for(int n = 0; n < 3; n++)
	{
		BackwardPlanU[n] = fftw_plan_dft_c2r_3d
						(Nx, Ny, Nz,
						 reinterpret_cast<fftw_complex*> (rcU[n]),
						 rlU[n],
						 FFTW_PATIENT);
	}
    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void ElasticitySolverSpectral::ReInitialize(ElasticProperties& EP)
{
    Nx  =  EP.Nx;
    Ny  =  EP.Ny;
    Nz  =  EP.Nz;
    Nz2 =     Nz/2+1;
    rlSIZE = Nx*Ny*Nz;
    rcSIZE = Nx*Ny*Nz2;

    dx = EP.dx;

    DPi_Nx = 2.0*Pi/double(Nx);
    DPi_Ny = 2.0*Pi/double(Ny);
    DPi_Nz = 2.0*Pi/double(Nz);

    Norm = 1.0/double(rlSIZE);
    // Arrays reallocation:
    for(int n = 0; n < 6; n++)
    {
        fftw_destroy_plan(ForwardPlanRHS[n]);

        delete[] rlRHSide[n];
        delete[] rcRHSide[n];
    }
    for(int n = 0; n < 3; n++)
    {
        fftw_destroy_plan(BackwardPlanU[n]);

    	delete[] rlU[n];
    	delete[] rcU[n];
        delete[] Q[n];
    }
    for(int n = 0; n < 9; n++)
    {
        fftw_destroy_plan(BackwardPlanDefGrad[n]);

        delete[] rlDefGrad[n];
        delete[] rcDefGrad[n];
    }

    for(int n = 0; n < 6; n++)
    {
        rlRHSide[n]  = new double[rlSIZE] ();
        memset(rlRHSide[n], 0, rlSIZE*sizeof(double));

        rcRHSide[n]  = new complex<double> [rcSIZE] ();
        memset(rcRHSide[n], 0, rcSIZE*sizeof(complex<double>));
    }
    for(int n = 0; n < 3; n++)
    {
    	rlU[n] = new double[rlSIZE] ();
        memset(rlU[n], 0, rlSIZE*sizeof(double));

		rcU[n] = new complex<double> [rcSIZE] ();
        memset(rcU[n], 0, rcSIZE*sizeof(complex<double>));

        Q[n]   = new double[rcSIZE] ();
        memset(Q[n], 0, rcSIZE*sizeof(double));
    }
    for(int n = 0; n < 9; n++)
    {
        rlDefGrad[n] = new double[rlSIZE] ();
        memset(rlDefGrad[n], 0, rlSIZE*sizeof(double));

        rcDefGrad[n] = new complex<double> [rcSIZE] ();
        memset(rcDefGrad[n], 0, rcSIZE*sizeof(complex<double>));
    }

    SetQ();

    for(int n = 0; n < 6; n++)
    {
        ForwardPlanRHS[n]  = fftw_plan_dft_r2c_3d
                        (Nx, Ny, Nz, rlRHSide[n],
                         reinterpret_cast<fftw_complex*> (rcRHSide[n]),
                         FFTW_MEASURE);
    }

    for(int n = 0; n < 9; n++)
    {
        BackwardPlanDefGrad[n] = fftw_plan_dft_c2r_3d
                             (Nx, Ny, Nz,
                              reinterpret_cast<fftw_complex*> (rcDefGrad[n]),
                              rlDefGrad[n],
                              FFTW_MEASURE);
    }

    for(int n = 0; n < 3; n++)
	{
		BackwardPlanU[n] = fftw_plan_dft_c2r_3d
						(Nx, Ny, Nz,
						 reinterpret_cast<fftw_complex*> (rcU[n]),
						 rlU[n],
						 FFTW_PATIENT);
	}
    Info::WriteStandard(thisclassname, "Reinitialized");
}

//++++++++++++ Destructor +++++++++++++++++++++++++

ElasticitySolverSpectral::~ElasticitySolverSpectral(void)
{
    for(int n = 0; n < 6; n++)
    {
        fftw_destroy_plan(ForwardPlanRHS[n]);

        delete[] rlRHSide[n];
        delete[] rcRHSide[n];
    }

    for(int n = 0; n < 3; n++)
    {
        fftw_destroy_plan(BackwardPlanU[n]);

    	delete[] rlU[n];
		delete[] rcU[n];
        delete[] Q[n];
    }
    for(int n = 0; n < 9; n++)
    {
        fftw_destroy_plan(BackwardPlanDefGrad[n]);
        delete[] rcDefGrad[n];
        delete[] rlDefGrad[n];
    }
//    fftw_cleanup_threads();
}

int ElasticitySolverSpectral::SetReferenceStrain(ElasticProperties& EP,
	BoundaryConditions& BC, double StrainAccuracy, double StressAccuracy,
	int MAXIterations, double dt)
{
	dMatrix6x6 Cij;                                                             ///< Starting elastic constants values
	Cij = (EP.MAXElasticConstants*1.1);

	vStrain  MAXStrainDeviation;
	vStress  AverageStress;
	vStress  TargetStress;

	MAXStrainDeviation.set_to_zero();
	AverageStress.set_to_zero();
	TargetStress.set_to_zero();

	vStrain   oldTargetStrain;
	oldTargetStrain.set_to_zero();
	EP.AverageStrain.set_to_zero();

	int    IterationCount = 0;
	double MAXStrainDifference = 0.0;
	double MAXTargetStrainDifference = 0.0;

	if (!EP.LargeDeformations)
	{
		for (int n = 0; n < 3; n++)
			EP.EffectiveAppliedStrain += EP.AppliedStrainRate*dt;
	}

	do // Iteration loop begin
	{
		/*OMP_PARALLEL_FOR_LOOP_BEGIN(XYZ, 0, rlSIZE, )
		{
			for (int n = 0; n < 9; n++)
			{
				rlDefGrad[n][XYZ] = 0.0;
			}
		}
		OMP_PARALLEL_FOR_LOOP_END*/

		IterationCount++;

		MAXStrainDifference = 0.0;
		MAXTargetStrainDifference = 0.0;
		MAXStrainDeviation.set_to_zero();
		AverageStress.set_to_zero();

		CalculateRHS(EP, Cij);
		ExecuteForwardFFT();
		CalculateFourierSolution(EP, Cij, false);
		ExecuteBackwardFFT(false);
		SetElasticProperties1(EP, MAXStrainDeviation, AverageStress);

		for (int n = 0; n < 6; n++)
		{
			if (MAXStrainDeviation[n] > MAXStrainDifference)
			{
				MAXStrainDifference = MAXStrainDeviation[n];
			}
			AverageStress[n] *= Norm;
		}

		TargetStress[0] = -AverageStress[0];
		TargetStress[1] = -AverageStress[1];
		TargetStress[2] = -AverageStress[2];

		oldTargetStrain = EP.AverageStrain;

		EP.AverageStrain.set_to_zero();

		for (int n = 0; n < 3; n++)
		for (int m = 0; m < 6; m++)
		{
			EP.AverageStrain[n] += TargetStress[m] * EP.AverageCompliences(n, m);
		}

//        for(int n = 0; n < 3; n++)
//        if(EP.AppStrainMask[n])
//        {
//            EP.AverageStrain[n] = EP.EffectiveAppliedStrain[n];
//        }

		SetElasticProperties2(EP, MAXStrainDeviation);

		for (int n = 0; n < 3; n++)
		if (fabs(EP.AverageStrain[n] - oldTargetStrain[n]) > MAXTargetStrainDifference)
		{
			MAXTargetStrainDifference = fabs(EP.AverageStrain[n] - oldTargetStrain[n]);
		}

		if (IterationCount > MAXIterations)
		{
			std::string message = "Maximum number of iterations (" + std::to_string(MAXIterations) + ") reached\n";
			message += Info::GetStandard("Strains converged to", Info::to_string_with_precision(MAXStrainDifference));
			Info::WriteWarning(message, thisclassname, "SetReferenceStrain()");
			break;
		}
	} // Iteration loop end
	while (MAXStrainDifference > StrainAccuracy or
		MAXTargetStrainDifference > StrainAccuracy);

	for (int n = 0; n < 3; n++)
	if (EP.AppStrainMask[n])
	{
		if (EP.LargeDeformations)
		{
			EP.AppliedStrain[n] += EP.AverageStrain[n];
		}
		else
		{
			EP.EffectiveAppliedStrain[n] += EP.AverageStrain[n];
		}
	}

	EP.EffectiveAppliedStrain += EP.AppliedStrainRate * dt;

	EP.SetBoundaryConditions(BC);

	return IterationCount;
}

int ElasticitySolverSpectral::Solve(ElasticProperties& EP,
         BoundaryConditions& BC, double StrainAccuracy, double StressAccuracy,
         int MAXIterations, double dt, bool getU)
{
    dMatrix6x6 Cij;                                                             ///< Starting elastic constants values
    Cij = (EP.MAXElasticConstants*1.1);

    vStrain  MAXStrainDeviation;
    vStress  AverageStress;
    vStress  TargetStress;

    MAXStrainDeviation.set_to_zero();
    AverageStress.set_to_zero();
    TargetStress.set_to_zero();

    vStrain   oldTargetStrain;
    oldTargetStrain.set_to_zero();
    EP.AverageStrain.set_to_zero();

    int    IterationCount = 0;
    double MAXStrainDifference = 0.0;
    double MAXTargetStrainDifference = 0.0;

	if (!EP.LargeDeformations)
	{
		for (int n = 0; n < 3; n++)
			EP.EffectiveAppliedStrain += EP.AppliedStrainRate*dt;
	}

    do // Iteration loop begin
    {
        /*OMP_PARALLEL_FOR_LOOP_BEGIN(XYZ,0,rlSIZE,)
        {
            for(int n = 0; n < 9; n++)
            {
                rlDefGrad[n][XYZ] = 0.0;
            }
        }
        OMP_PARALLEL_FOR_LOOP_END*/

        IterationCount++;

        MAXStrainDifference = 0.0;
        MAXTargetStrainDifference = 0.0;
        MAXStrainDeviation.set_to_zero();
        AverageStress.set_to_zero();

        CalculateRHS(EP, Cij);
        ExecuteForwardFFT();
        CalculateFourierSolution(EP, Cij, getU);
        ExecuteBackwardFFT(getU);
        SetElasticProperties1(EP, MAXStrainDeviation, AverageStress);

        for(int n = 0; n < 6; n++)
        {
            if(MAXStrainDeviation[n] > MAXStrainDifference)
            {
               MAXStrainDifference = MAXStrainDeviation[n];
            }
            AverageStress[n] *= Norm;
        }

        TargetStress[0] = -AverageStress[0] + EP.EffectiveAppliedStress[0];
        TargetStress[1] = -AverageStress[1] + EP.EffectiveAppliedStress[1];
        TargetStress[2] = -AverageStress[2] + EP.EffectiveAppliedStress[2];

        oldTargetStrain = EP.AverageStrain;

        EP.AverageStrain.set_to_zero();

        for(int n = 0; n < 3; n++)
        for(int m = 0; m < 6; m++)
        {
            EP.AverageStrain[n] += TargetStress[m]*EP.AverageCompliences(n,m);
        }
        if (EP.KeepAspectRatio)
		{
			double trace = (1.0/3.0)*(EP.AverageStrain[0] + EP.AverageStrain[1] + EP.AverageStrain[2]);
			EP.AverageStrain[0] = trace;
			EP.AverageStrain[1] = trace;
			EP.AverageStrain[2] = trace;
		}

        SetElasticBoundaryConditions(EP, TargetStress);

//        for(int n = 0; n < 3; n++)
//        if(EP.AppStrainMask[n])
//        {
//            EP.AverageStrain[n] = EP.EffectiveAppliedStrain[n];
//        }

        SetElasticProperties2(EP, MAXStrainDeviation);

        for(int n = 0; n < 3; n++)
        if(fabs(EP.AverageStrain[n] - oldTargetStrain[n]) > MAXTargetStrainDifference)
        {
            MAXTargetStrainDifference = fabs(EP.AverageStrain[n] - oldTargetStrain[n]);
        }

        if(IterationCount > MAXIterations)
        {
            std::string message = "Maximum number of iterations (" + std::to_string(MAXIterations) + ") reached\n";
               message += Info::GetStandard("Strains converged to", Info::to_string_with_precision(MAXStrainDifference));
               Info::WriteWarning(message, thisclassname, "Solve()");
               break;
        }
    } // Iteration loop end
    while(MAXStrainDifference > StrainAccuracy or
          MAXTargetStrainDifference > StrainAccuracy);

	EP.EffectiveAppliedStrain += EP.AppliedStrainRate * dt;

    EP.SetBoundaryConditions(BC);

    return IterationCount;
}

int ElasticitySolverSpectral::Solve(ElasticProperties& EP, Orientations& OR,
         BoundaryConditions& BC, double StrainAccuracy, double StressAccuracy,
         int MAXIterations, double dt, bool getU)
{
    dMatrix6x6 Cij;                                                             ///< Starting elastic constants values
    Cij = (EP.MAXElasticConstants*1.1);

    vStrain MAXStrainDeviation;
    vStress AverageStress;
    vStress TargetStress;

    AverageStress.set_to_zero();
    MAXStrainDeviation.set_to_zero();
    TargetStress.set_to_zero();

    vStrain   oldTargetStrain;
    oldTargetStrain.set_to_zero();
    EP.AverageStrain.set_to_zero();

    int    IterationCount = 0;
    double MAXStrainDifference = 0.0;
    double MAXTargetStrainDifference = 0.0;

	if (!EP.LargeDeformations)
	{
		for (int n = 0; n < 3; n++)
			EP.EffectiveAppliedStrain += EP.AppliedStrainRate*dt;
	}

    do // Iteration loop begin
    {
        /*OMP_PARALLEL_FOR_LOOP_BEGIN(XYZ,0,rlSIZE,)
        {
            for(int n = 0; n < 9; n++)
            {
                rlDefGrad[n][XYZ] = 0.0;
            }
        }
        OMP_PARALLEL_FOR_LOOP_END*/

        IterationCount++;

        MAXStrainDifference = 0.0;
        MAXTargetStrainDifference = 0.0;
        MAXStrainDeviation.set_to_zero();
        AverageStress.set_to_zero();

        CalculateRHS(EP, Cij);
        ExecuteForwardFFT();
        CalculateFourierSolution(EP, Cij, getU);
        ExecuteBackwardFFT(getU);
        SetElasticProperties1(EP, MAXStrainDeviation, AverageStress);

        for(int n = 0; n < 6; n++)
        {
            if(MAXStrainDeviation[n] > MAXStrainDifference)
            {
               MAXStrainDifference = MAXStrainDeviation[n];
            }
            AverageStress[n] *= Norm;
        }

        TargetStress[0] = -AverageStress[0] + EP.EffectiveAppliedStress[0];
        TargetStress[1] = -AverageStress[1] + EP.EffectiveAppliedStress[1];
        TargetStress[2] = -AverageStress[2] + EP.EffectiveAppliedStress[2];

        oldTargetStrain = EP.AverageStrain;

        EP.AverageStrain.set_to_zero();

        for(int n = 0; n < 3; n++)
        for(int m = 0; m < 6; m++)
        {
            EP.AverageStrain[n] += TargetStress[m]*EP.AverageCompliences(n,m);
        }

        if (EP.KeepAspectRatio)
        {
            double trace = (1.0/3.0)*(EP.AverageStrain[0] + EP.AverageStrain[1] + EP.AverageStrain[2]);
            EP.AverageStrain[0] = trace;
            EP.AverageStrain[1] = trace;
            EP.AverageStrain[2] = trace;
        }

        SetElasticBoundaryConditions(EP, TargetStress);
//        for(int n = 0; n < 3; n++)
//        if(EP.AppStrainMask[n])
//        {
//            EP.AverageStrain[n] = EP.EffectiveAppliedStrain[n];
//        }

        SetElasticProperties2(EP, MAXStrainDeviation);

        for(int n = 0; n < 3; n++)
        if(fabs(EP.AverageStrain[n] - oldTargetStrain[n]) > MAXTargetStrainDifference)
        {
            MAXTargetStrainDifference = fabs(EP.AverageStrain[n] - oldTargetStrain[n]);
        }

        if(IterationCount > MAXIterations)
        {
            std::string message = "Maximum number of iterations (" + std::to_string(MAXIterations) + ") reached\n";
               message += Info::GetStandard("Strains converged to", Info::to_string_with_precision(MAXStrainDifference));
               Info::WriteWarning(message, thisclassname, "Solve()");
               break;
        }
    } // Iteration loop end
    while(MAXStrainDifference > StrainAccuracy or
          MAXTargetStrainDifference > StrainAccuracy);

    if(EP.LargeDeformations)
    {
        EP.SetVelocityGradient(rlDefGrad, dt);
        EP.SetRotations(OR, dt);
    }

    EP.SetBoundaryConditions(BC);
    OR.SetBoundaryConditions(BC);

    return IterationCount;
}

void ElasticitySolverSpectral::CalculateRHS(ElasticProperties& EP, dMatrix6x6& Cij)
{
	dMatrix3x3 unity;
	unity.set_to_unity();

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Strains, 0,)
	{
        dMatrix3x3 locDefGrad;
        locDefGrad(0,0) = rlDefGrad[0][k + Nz*(j + Ny*i)] + 1.0;
        locDefGrad(1,1) = rlDefGrad[4][k + Nz*(j + Ny*i)] + 1.0;
        locDefGrad(2,2) = rlDefGrad[8][k + Nz*(j + Ny*i)] + 1.0;

        locDefGrad(0,1) = rlDefGrad[1][k + Nz*(j + Ny*i)];
        locDefGrad(0,2) = rlDefGrad[2][k + Nz*(j + Ny*i)];
        locDefGrad(1,2) = rlDefGrad[5][k + Nz*(j + Ny*i)];

        locDefGrad(1,0) = rlDefGrad[3][k + Nz*(j + Ny*i)];
        locDefGrad(2,0) = rlDefGrad[6][k + Nz*(j + Ny*i)];
        locDefGrad(2,1) = rlDefGrad[7][k + Nz*(j + Ny*i)];

        dMatrix3x3 locStrainTensorFinite = (locDefGrad.transposed()*locDefGrad - unity)*0.5;
        dMatrix3x3 locStrainTensorSmall = (locDefGrad.transposed() + locDefGrad)*0.5 - unity;
        vStrain locStrainDiff = locStrainTensorFinite.VoigtStrain() - locStrainTensorSmall.VoigtStrain();
		for(int n = 0; n <  6; n++)
		{
			rlRHSide[n][k + Nz*(j + Ny*i)] = 0.0;
			for(int m = 0; m < 6; m++)
			{
				rlRHSide[n][k + Nz*(j + Ny*i)] +=
				EP.EffectiveElasticConstants(i, j, k)(n,m)*
				(EP.EffectiveEigenStrains(i, j, k)[m] - locStrainDiff[m])-
				(EP.EffectiveElasticConstants(i, j, k)(n,m) - Cij(n,m))*(EP.StrainIncrements(i,j,k)[m]);
			}
		}
	}
    OMP_PARALLEL_STORAGE_LOOP_END
    
    if(Nx > EP.Nx)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Strains, 0,)
        for(int n = 0; n < 6; n++)
        {
            int x = Nx - i - 1;
            for(int m = 0; m < 6; m++)
            {
                rlRHSide[n][k + Nz*(j + Ny*x)] = rlRHSide[n][k + Nz*(j + Ny*i)];
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    if(Ny > EP.Ny)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Strains, 0,)
        for(int n = 0; n < 6; n++)
        {
            int y = Ny - j - 1;
            for(int m = 0; m < 6; m++)
            {
                rlRHSide[n][k + Nz*(y + Ny*i)] = rlRHSide[n][k + Nz*(j + Ny*i)];
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    if(Nz > EP.Nz)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Strains, 0,)
        for(int n = 0; n < 6; n++)
        {
            int z = Nz - k - 1;
            for(int m = 0; m < 6; m++)
            {
                rlRHSide[n][z + Nz*(j + Ny*i)] = rlRHSide[n][k + Nz*(j + Ny*i)];
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

void ElasticitySolverSpectral::CalculateFourierSolution(ElasticProperties& EP, dMatrix6x6 Cij, bool getU)
{
    double Norm = 1.0/double(rlSIZE);

    OMP_PARALLEL_FOR_LOOP_BEGIN(XYZ,0,rcSIZE,)
    {
        double Qx = Q[0][XYZ];
        double Qy = Q[1][XYZ];
        double Qz = Q[2][XYZ];

        complex<double> rhsX = -I*(Qx*rcRHSide[0][XYZ] +
                                   Qy*rcRHSide[5][XYZ] +
                                   Qz*rcRHSide[4][XYZ]);
        complex<double> rhsY = -I*(Qx*rcRHSide[5][XYZ] +
                                   Qy*rcRHSide[1][XYZ] +
                                   Qz*rcRHSide[3][XYZ]);
        complex<double> rhsZ = -I*(Qx*rcRHSide[4][XYZ] +
                                   Qy*rcRHSide[3][XYZ] +
                                   Qz*rcRHSide[2][XYZ]);

        double a11 = (Cij(0,0)*Qx*Qx + 2.0*Cij(0,5)*Qx*Qy + Cij(5,5)*Qy*Qy +
                  2.0*Cij(0,4)*Qx*Qz + 2.0*Cij(4,5)*Qy*Qz + Cij(4,4)*Qz*Qz);

        double a21 = (Cij(0,5)*Qx*Qx + Cij(0,1)*Qx*Qy + Cij(5,5)*Qx*Qy +
                      Cij(1,5)*Qy*Qy + Cij(0,3)*Qx*Qz + Cij(4,5)*Qx*Qz +
                      Cij(1,4)*Qy*Qz + Cij(3,5)*Qy*Qz + Cij(3,4)*Qz*Qz);

        double a31 = (Cij(0,4)*Qx*Qx + Cij(0,3)*Qx*Qy + Cij(4,5)*Qx*Qy +
                      Cij(3,5)*Qy*Qy + Cij(0,2)*Qx*Qz + Cij(4,4)*Qx*Qz +
                      Cij(2,5)*Qy*Qz + Cij(3,4)*Qy*Qz + Cij(2,4)*Qz*Qz);

        double a12 = a21;

        double a22 = (Cij(5,5)*Qx*Qx + 2.0*Cij(1,5)*Qx*Qy + Cij(1,1)*Qy*Qy +
                  2.0*Cij(3,5)*Qx*Qz + 2.0*Cij(1,3)*Qy*Qz + Cij(3,3)*Qz*Qz);

        double a32 = (Cij(4,5)*Qx*Qx + Cij(1,4)*Qx*Qy + Cij(3,5)*Qx*Qy +
                      Cij(1,3)*Qy*Qy + Cij(2,5)*Qx*Qz + Cij(3,4)*Qx*Qz +
                      Cij(1,2)*Qy*Qz + Cij(3,3)*Qy*Qz + Cij(2,3)*Qz*Qz);

        double a13 = a31;

        double a23 = a32;

        double a33 = (Cij(4,4)*Qx*Qx + 2.0*Cij(3,4)*Qx*Qy + Cij(3,3)*Qy*Qy +
                  2.0*Cij(2,4)*Qx*Qz + 2.0*Cij(2,3)*Qy*Qz + Cij(2,2)*Qz*Qz);

        double denominator = (-a13*a22*a31 + a12*a23*a31 + a13*a21*a32 -
                               a11*a23*a32 - a12*a21*a33 + a11*a22*a33);

        if(std::abs(denominator) > DBL_EPSILON and std::abs(denominator) < DBL_MAX)
        {
            denominator = 1.0/denominator;
        }
        else
        {
            denominator = 0.0;
        }

        complex<double> locUrcX = (-a23*a32*rhsX + a22*a33*rhsX + a13*a32*rhsY -
                                    a12*a33*rhsY - a13*a22*rhsZ + a12*a23*rhsZ)*denominator*Norm;

        complex<double> locUrcY = ( a23*a31*rhsX - a21*a33*rhsX - a13*a31*rhsY +
                                    a11*a33*rhsY + a13*a21*rhsZ - a11*a23*rhsZ)*denominator*Norm;

        complex<double> locUrcZ = (-a22*a31*rhsX + a21*a32*rhsX + a12*a31*rhsY -
                                    a11*a32*rhsY - a12*a21*rhsZ + a11*a22*rhsZ)*denominator*Norm;

        if(getU)
        {
        	rcU[0][XYZ] = locUrcX;
        	rcU[1][XYZ] = locUrcY;
        	rcU[2][XYZ] = locUrcZ;
        }
        rcDefGrad[0][XYZ] = I*(Qx*locUrcX);
        rcDefGrad[1][XYZ] = I*(Qy*locUrcX);
        rcDefGrad[2][XYZ] = I*(Qz*locUrcX);
        rcDefGrad[3][XYZ] = I*(Qx*locUrcY);
        rcDefGrad[4][XYZ] = I*(Qy*locUrcY);
        rcDefGrad[5][XYZ] = I*(Qz*locUrcY);
        rcDefGrad[6][XYZ] = I*(Qx*locUrcZ);
        rcDefGrad[7][XYZ] = I*(Qy*locUrcZ);
        rcDefGrad[8][XYZ] = I*(Qz*locUrcZ);
    }
    OMP_PARALLEL_FOR_LOOP_END
}

void ElasticitySolverSpectral::SetElasticProperties1(ElasticProperties& EP,
                    vStrain& MAXStrainDeviation, vStress& AverageStress)
{
	dMatrix3x3 unity;
	unity.set_to_unity();
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.StrainIncrements,0,/*reduction(dVector6MaxComp:MAXStrainDeviation)*/)
    {
        dMatrix3x3 locDefGrad;
        locDefGrad(0,0) = rlDefGrad[0][k + Nz*(j + Ny*i)] + 1.0;
        locDefGrad(1,1) = rlDefGrad[4][k + Nz*(j + Ny*i)] + 1.0;
        locDefGrad(2,2) = rlDefGrad[8][k + Nz*(j + Ny*i)] + 1.0;

        locDefGrad(0,1) = rlDefGrad[1][k + Nz*(j + Ny*i)];
        locDefGrad(0,2) = rlDefGrad[2][k + Nz*(j + Ny*i)];
        locDefGrad(1,2) = rlDefGrad[5][k + Nz*(j + Ny*i)];

        locDefGrad(1,0) = rlDefGrad[3][k + Nz*(j + Ny*i)];
        locDefGrad(2,0) = rlDefGrad[6][k + Nz*(j + Ny*i)];
        locDefGrad(2,1) = rlDefGrad[7][k + Nz*(j + Ny*i)];

        dMatrix3x3 locStrainTensor = (locDefGrad.transposed()*locDefGrad - unity)*0.5;
        //dMatrix3x3 locStrainTensor = (locDefGrad.transposed() + locDefGrad)*0.5 - unity;
        vStrain locStrain = locStrainTensor.VoigtStrain();

        for(int n = 0; n < 6; n++)
        {
            double locStrainDifference = fabs(EP.StrainIncrements(i,j,k)[n] - locStrain[n] - EP.AverageStrain[n]);
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                if(locStrainDifference > MAXStrainDeviation[n])
                {
                    MAXStrainDeviation[n] = locStrainDifference;
                }
            }
        }
        for(int n = 0; n < 6; n++)
        {
            double locStress = 0.0;
            for(int m = 0; m < 6; m++)
            {
                locStress += EP.EffectiveElasticConstants(i, j, k)(n,m)*
                          (locStrain[m] - EP.EffectiveEigenStrains(i, j, k)[m]);// + EP.RemeshedStrain[m]);
            }
            #ifdef _OPENMP
            #pragma omp atomic
            #endif
            AverageStress[n] += locStress;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySolverSpectral::SetElasticProperties2(ElasticProperties& EP,
                                                   vStrain& MAXStrainDeviation)
{
	dMatrix3x3 unity;
	unity.set_to_unity();
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.StrainIncrements,0,/*reduction(dVector6MaxComp:MAXStrainDeviation)*/)
    {
		dMatrix3x3 locDefGrad;
        locDefGrad(0,0) = rlDefGrad[0][k + Nz*(j + Ny*i)] + 1.0;
        locDefGrad(1,1) = rlDefGrad[4][k + Nz*(j + Ny*i)] + 1.0;
        locDefGrad(2,2) = rlDefGrad[8][k + Nz*(j + Ny*i)] + 1.0;

        locDefGrad(0,1) = rlDefGrad[1][k + Nz*(j + Ny*i)];
        locDefGrad(0,2) = rlDefGrad[2][k + Nz*(j + Ny*i)];
        locDefGrad(1,2) = rlDefGrad[5][k + Nz*(j + Ny*i)];

        locDefGrad(1,0) = rlDefGrad[3][k + Nz*(j + Ny*i)];
        locDefGrad(2,0) = rlDefGrad[6][k + Nz*(j + Ny*i)];
        locDefGrad(2,1) = rlDefGrad[7][k + Nz*(j + Ny*i)];

		dMatrix3x3 locStrainTensor = (locDefGrad.transposed()*locDefGrad - unity)*0.5;
        //dMatrix3x3 locStrainTensor = (locDefGrad.transposed() + locDefGrad)*0.5 - unity;
        vStrain locStrain = locStrainTensor.VoigtStrain();

		for(int n = 0; n < 6; n++)
		{
			locStrain[n] += EP.AverageStrain[n];

            double locStrainDifference = fabs(EP.StrainIncrements(i,j,k)[n] - locStrain[n]);

            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                if(locStrainDifference > MAXStrainDeviation[n])
                {
                    MAXStrainDeviation[n] = locStrainDifference;
                }
            }
            EP.StrainIncrements(i,j,k)[n] = locStrain[n];
            if(!EP.LargeDeformations) EP.Strains(i,j,k)[n] = locStrain[n];
        }
        for(int n = 0; n < 6; n++)
        {
            double locStress = 0.0;
            for(int m = 0; m < 6; m++)
            {
                locStress += EP.EffectiveElasticConstants(i, j, k)(n,m)*
                          (locStrain[m] - EP.EffectiveEigenStrains(i, j, k)[m]);// + EP.RemeshedStrain[m]);
            }
            if(EP.LargeDeformations)
            {
                EP.StressIncrements(i,j,k)[n] = locStress;
            }
            else
            {
                EP.Stresses(i,j,k)[n] = locStress;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticitySolverSpectral::ExecuteForwardFFT(void)
{
    #pragma omp parallel sections// OMP BEGIN
    {
        #pragma omp section
        {
            fftw_execute(ForwardPlanRHS[0]);
        }
        #pragma omp section
        {
            fftw_execute(ForwardPlanRHS[1]);
        }
        #pragma omp section
        {
            fftw_execute(ForwardPlanRHS[2]);
        }
        #pragma omp section
        {
            fftw_execute(ForwardPlanRHS[3]);
        }
        #pragma omp section
        {
            fftw_execute(ForwardPlanRHS[4]);
        }
        #pragma omp section
        {
            fftw_execute(ForwardPlanRHS[5]);
        }
    }//OMP END
}

void ElasticitySolverSpectral::ExecuteBackwardFFT(bool getU)
{
    #pragma omp parallel sections // OMP BEGIN
    {
        #pragma omp section
        {
            fftw_execute(BackwardPlanDefGrad[0]);
        }
        #pragma omp section
        {
            fftw_execute(BackwardPlanDefGrad[1]);
        }
        #pragma omp section
        {
            fftw_execute(BackwardPlanDefGrad[2]);
        }
        #pragma omp section
        {
            fftw_execute(BackwardPlanDefGrad[3]);
        }
        #pragma omp section
        {
            fftw_execute(BackwardPlanDefGrad[4]);
        }
        #pragma omp section
        {
            fftw_execute(BackwardPlanDefGrad[5]);
        }
        #pragma omp section
        {
            fftw_execute(BackwardPlanDefGrad[6]);
        }
        #pragma omp section
        {
            fftw_execute(BackwardPlanDefGrad[7]);
        }
        #pragma omp section
        {
            fftw_execute(BackwardPlanDefGrad[8]);
        }
        #pragma omp section
		{
			if(getU) fftw_execute(BackwardPlanU[0]);
		}
		#pragma omp section
		{
			if(getU) fftw_execute(BackwardPlanU[1]);
		}
		#pragma omp section
		{
			if(getU) fftw_execute(BackwardPlanU[2]);
		}
    }
}

void ElasticitySolverSpectral::SetQ()
{
    for(int i = 0; i < Nx ; i++)
    for(int j = 0; j < Ny ; j++)
    for(int k = 0; k < Nz2; k++)
    {
         int XYZ = k + Nz2*(j + Ny*i);

         Q[0][XYZ] = DPi_Nx*(i*(i <= Nx/2) - (Nx-i)*(i > Nx/2))/dx;
         Q[1][XYZ] = DPi_Ny*(j*(j <= Ny/2) - (Ny-j)*(j > Ny/2))/dx;
         Q[2][XYZ] = DPi_Nz*(k*(k <= Nz/2) - (Nz-k)*(k > Nz/2))/dx;
    }
}

void ElasticitySolverSpectral::SetElasticBoundaryConditions(ElasticProperties& EP, vStress TargetStress)
{
    if (EP.AverageElasticConstants(0,0) == 0 or EP.AverageElasticConstants(1,1) == 0 or EP.AverageElasticConstants(2,2) == 0
            or EP.AverageElasticConstants(1,2) == 0  or EP.AverageElasticConstants(0,2) == 0  or EP.AverageElasticConstants(0,1) == 0)
    {
        std::string message = "Zero component in EP.AverageElasticConstants.";
           Info::WriteWarning(message, thisclassname, "SetElasticBoundaryConditions()");
           exit(1);
    }

    if(EP.AppStrainMask[0] and !EP.AppStrainMask[1] and !EP.AppStrainMask[2])   /// XX
    {
		EP.AverageStrain[0] = EP.EffectiveAppliedStrain[0];// +EP.AppliedStrainRate[0] * EP.dt;
        EP.AverageStrain[1] = (- EP.AverageElasticConstants(0,2) * EP.AverageElasticConstants(1,2) * EP.EffectiveAppliedStrain[0] +
                                EP.AverageElasticConstants(0,1) * EP.AverageElasticConstants(2,2) * EP.EffectiveAppliedStrain[0] -
                                EP.AverageElasticConstants(2,2) * TargetStress[1] + EP.AverageElasticConstants(1,2) * TargetStress[2]) /
                                (EP.AverageElasticConstants(1,2) * EP.AverageElasticConstants(1,2) - EP.AverageElasticConstants(1,1) * EP.AverageElasticConstants(2,2));
//        EP.AverageStrain[2] = (- EP.AverageElasticConstants(0,1) * EP.AverageStrain[0] - EP.AverageElasticConstants(1,1) * EP.EffectiveAppliedStrain[1] + TargetStress[1]) /
//                                EP.AverageElasticConstants(1,2);
        EP.AverageStrain[2] = (- EP.AverageElasticConstants(1,1) * EP.AverageStrain[1] - EP.AverageElasticConstants(0,1) * EP.EffectiveAppliedStrain[0] + TargetStress[1]) /
                                EP.AverageElasticConstants(1,2);
    }

    if(!EP.AppStrainMask[0] and EP.AppStrainMask[1] and !EP.AppStrainMask[2])   /// YY
    {
        EP.AverageStrain[0] = (- EP.AverageElasticConstants(0,2) * EP.AverageElasticConstants(1,2) * EP.EffectiveAppliedStrain[1] +
                                EP.AverageElasticConstants(0,1) * EP.AverageElasticConstants(2,2) * EP.EffectiveAppliedStrain[1] -
                                EP.AverageElasticConstants(2,2) * TargetStress[0] + EP.AverageElasticConstants(0,2) * TargetStress[2]) /
                                (EP.AverageElasticConstants(0,2) * EP.AverageElasticConstants(0,2) - EP.AverageElasticConstants(0,0) * EP.AverageElasticConstants(2,2));
        EP.AverageStrain[1] = EP.EffectiveAppliedStrain[1];// + EP.AppliedStrainRate[1] * EP.dt;
        EP.AverageStrain[2] = (- EP.AverageElasticConstants(0,0) * EP.AverageStrain[0] - EP.AverageElasticConstants(0,1) * EP.EffectiveAppliedStrain[1] + TargetStress[0]) /
                                EP.AverageElasticConstants(0,2);
    }

    if(!EP.AppStrainMask[0] and !EP.AppStrainMask[1] and EP.AppStrainMask[2])   /// ZZ
    {
        EP.AverageStrain[0] = (EP.AverageElasticConstants(0,2) * EP.AverageElasticConstants(1,1) * EP.EffectiveAppliedStrain[2] -
                                EP.AverageElasticConstants(0,1) * EP.AverageElasticConstants(1,2) * EP.EffectiveAppliedStrain[2] -
                                EP.AverageElasticConstants(1,1) * TargetStress[0] + EP.AverageElasticConstants(0,1) * TargetStress[1]) /
                                (EP.AverageElasticConstants(0,1) * EP.AverageElasticConstants(0,1) - EP.AverageElasticConstants(0,0) * EP.AverageElasticConstants(1,1));

        EP.AverageStrain[1] = (- EP.AverageElasticConstants(0,0) * EP.AverageStrain[0] - EP.AverageElasticConstants(0,2) * EP.EffectiveAppliedStrain[2] + TargetStress[0]) /
                                EP.AverageElasticConstants(0,1);
		EP.AverageStrain[2] = EP.EffectiveAppliedStrain[2];// +EP.AppliedStrainRate[2] * EP.dt;
    }

    if(EP.AppStrainMask[0] and EP.AppStrainMask[1] and !EP.AppStrainMask[2])    /// XX & YY
    {
        EP.AverageStrain[0] = EP.EffectiveAppliedStrain[0];// + EP.AppliedStrainRate[0] * EP.dt;
        EP.AverageStrain[1] = EP.EffectiveAppliedStrain[1];// + EP.AppliedStrainRate[1] * EP.dt;
        EP.AverageStrain[2] = (- EP.AverageElasticConstants(0,2) * EP.EffectiveAppliedStrain[0] - EP.AverageElasticConstants(1,2) * EP.EffectiveAppliedStrain[1] + TargetStress[2]) /
                                EP.AverageElasticConstants(2,2);
    }

    if(EP.AppStrainMask[0] and !EP.AppStrainMask[1] and EP.AppStrainMask[2])    /// XX & ZZ
    {
        EP.AverageStrain[0] = EP.EffectiveAppliedStrain[0];// + EP.AppliedStrainRate[0] * EP.dt;
        EP.AverageStrain[1] = (- EP.AverageElasticConstants(0,1) * EP.EffectiveAppliedStrain[0] - EP.AverageElasticConstants(1,2) * EP.EffectiveAppliedStrain[2] + TargetStress[1]) /
                                EP.AverageElasticConstants(1,1);
        EP.AverageStrain[2] = EP.EffectiveAppliedStrain[2];// + EP.AppliedStrainRate[2] * EP.dt;
    }

    if(!EP.AppStrainMask[0] and EP.AppStrainMask[1] and EP.AppStrainMask[2])    /// YY & ZZ
    {
        EP.AverageStrain[0] = (- EP.AverageElasticConstants(0,1) * EP.EffectiveAppliedStrain[1] - EP.AverageElasticConstants(0,2) * EP.EffectiveAppliedStrain[2] + TargetStress[0]) /
                                EP.AverageElasticConstants(0,0);
        EP.AverageStrain[1] = EP.EffectiveAppliedStrain[1];// + EP.AppliedStrainRate[1] * EP.dt;
        EP.AverageStrain[2] = EP.EffectiveAppliedStrain[2];// + EP.AppliedStrainRate[2] * EP.dt;
    }

    if(EP.AppStrainMask[0] and EP.AppStrainMask[1] and EP.AppStrainMask[2])    /// YY & ZZ
    {
        EP.AverageStrain[0] = EP.EffectiveAppliedStrain[0];// + EP.AppliedStrainRate[0] * EP.dt;
        EP.AverageStrain[1] = EP.EffectiveAppliedStrain[1];// + EP.AppliedStrainRate[1] * EP.dt;
        EP.AverageStrain[2] = EP.EffectiveAppliedStrain[2];// + EP.AppliedStrainRate[2] * EP.dt;
    }
}
} // namespace openphase
