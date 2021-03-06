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

#ifndef ELASTICITYSOLVERSPECTRAL_H
#define ELASTICITYSOLVERSPECTRAL_H

/******************************************************************************
*  This module is based on the iterative algorithm of [S.Y.Hu, L.Q.Chen,      *
*  Acta. Mater. 49 (2001) 1879 - 1890] with modifications allowing to use     *
*  effective space dependent elastic constants not bound to concentration.    *
*  The boundary condition allowing free volume expansion is also introduced   *
*  following the approach outlined in [I.Steinbach, M.Apel, Physica D 217     *
*  (2006) 153 - 160].                                                         *
******************************************************************************/
#include "fftw3.h"
#include "Tools/Includes.h"
#include "Mechanics/ElasticitySolvers/ElasticitySolver.h"

namespace opensim
{
class ElasticProperties;
class Orientations;
class Settings;

class ElasticitySolverSpectral : public ElasticitySolver                        ///< Elastic problem solver based on spectral method.
{
 public:
    ElasticitySolverSpectral(){};
    ElasticitySolverSpectral(Settings& locSettings);
    ~ElasticitySolverSpectral(void);                                            ///< Destructor
    void Initialize(Settings& locSettings);                                     ///< Constructor
    void Initialize(Settings& locSettings, BoundaryConditions& BC);             ///< Constructor (considering non periodic boundary conditions)
    void ReInitialize(ElasticProperties& EP);                                   ///< Needed if the system has been remeshed

	int SetReferenceStrain(ElasticProperties& EP,
		BoundaryConditions& BC, double StrainAccuracy, double StressAccuracy,
		int MAXIterations, double dt);											///< Relaxes the system and sets (Effective)AppliedStrain to averagestrain
    int Solve(ElasticProperties& EP, Orientations& OR, BoundaryConditions& BC,
              double StrainAccuracy, double StressAccuracy, int MAXIterations,
              double dt, bool getU = false);                                    ///< Solves elastic problem
    int Solve(ElasticProperties& EP, BoundaryConditions& BC,
              double StrainAccuracy, double StressAccuracy, int MAXIterations,
              double dt, bool getU = false);                                    ///< Solves elastic problem
    
    void CalculateRHS(ElasticProperties& EP, dMatrix6x6& Cij);
    void ExecuteForwardFFT();
    void CalculateFourierSolution(ElasticProperties& EP, dMatrix6x6 Cij, bool getU);
    void ExecuteBackwardFFT(bool getU);
    void SetElasticProperties1(ElasticProperties& EP, vStrain& MAXStrainDeviation, 
                                                      vStress& AverageStress);
    void SetElasticProperties2(ElasticProperties& EP, vStrain& MAXStrainDeviation);
    void SetQ(void);                                                            ///< Sets the wave vectors

    void SetElasticBoundaryConditions(ElasticProperties& EP, vStress TargetStress);

    long int Nx;                                                                ///< System size along X direction
    long int Ny;                                                                ///< System size along y direction
    long int Nz;                                                                ///< System size along z direction

    double               * rlU[3];												///< Displacements in real space

 private:
    double dx;                                                                  ///< Grid spacing
    long int Nz2;                                                               ///< Half of the system size along Z direction

    long int rlSIZE;                                                            ///< System size (real space): Nx*Ny*Nz
    long int rcSIZE;                                                            ///< System size (Fourier space): Nx*Ny*(Nz/2+1)

    double Norm;                                                                ///< Normalization coefficient: 1.0/SIZE
    double DPi_Nx;                                                              ///< 2.0*Pi/Nx constant
    double DPi_Ny;                                                              ///< 2.0*Pi/Ny constant
    double DPi_Nz;                                                              ///< 2.0*Pi/Nz constant

    dMatrix6x6 C0;
    dMatrix6x6 C0inverse;

    double               * rlRHSide[6];                                         ///< Right hand side in real space
    std::complex<double> * rcRHSide[6];                                         ///< Right hand side in Fourier space

    std::complex<double> * rcU[3];												///< Displacements in reciprocal space

    double               * Q[3];                                                ///< Wave vectors

    double               * rlDefGrad[9];                                        ///< Deformation gradient entries in real space
    std::complex<double> * rcDefGrad[9];                                        ///< Deformation gradient entries in Fourier space

    fftw_plan  ForwardPlanRHS[6];                                               ///< Forward FFT plans for RHSide
    fftw_plan  BackwardPlanDefGrad[9];                                          ///< Backward FFT plans for deformation gradients
    fftw_plan  BackwardPlanU[3];                                                ///< Backward FFT plans for displacements
};
} // namespace openphase
#endif //SpectralElasticitySolver
