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

#ifndef NOISE_H
#define NOISE_H

#include "Tools/Includes.h"
#include "Tools/Node.h"
#include "fftw3.h"

namespace opensim
{
class DrivingForce;
class PhaseField;
class Settings;
class Temperature;

class Noise : public OPObject                                                   ///< Calculates a noise, which can be added to the DrivingForce
{
    public:
        Noise(void){};                                                          ///< Simple constructor, but Initialize and ReadInput have to be called additionally
        Noise(Settings& locSettings);                                           ///< Constructor includes the call of Initialize, bur ReadInput has to be called separately
        Noise(Settings& locSettings, std::string InputFileName);                ///< Constructor includes the call of Initialize and ReadInput

        void AddToDrivingForce(int tStep, Temperature& Temp, DrivingForce& DF); ///< Adds Fields to driving force
        void Initialize(Settings& locSettings);                                 ///< Allocates memory, initializes the settings
        void Initialize(Settings& locSettings, std::string InputFileName);      ///< Allocates memory, initializes the settings, reads input
        void ReadInput(std::string InputFileName);                              ///< Reads input from input file
        void WriteVTK(int tStep);                                               ///< Writes the noise acting on the driving force from phase field with indexB onto phase field with indexA

        double gamma;                                                           ///< Macroscopic response of the system
        int    BCells;                                                          ///< Number of boundary cells around the computation domain
        int    KCutOff;                                                         ///< Number of Fourier coefficient in each spatial direction that will be set to zero
        int    Nphases;                                                         ///< Number of phase-fields
        int    Nx;                                                              ///< Grid size in x-direction
        int    Ny;                                                              ///< Grid size in y-direction
        int    Nz;                                                              ///< Grid size in z-direction
        int    RandomSeed;                                                      ///< Defines the random number table (used to make results reproducible)
        int    TimeSteps;                                                       ///< Number of iteration after which a new noise will be scrabbled
        int    dt;                                                              ///< Size time discretization
        int    dx;                                                              ///< Grid spacing in x-direction
        int    dy;                                                              ///< Grid spacing in y-direction
        int    dz;                                                              ///< Grid spacing in z-direction
        //const double kBoltzmann;                                 ///< Physical constant

    protected:

        void Generate(void);                                                    ///< Generates a new noise field (it will be interpolated between new and old noise field)

        Storage3D< double, 0 >           Raw;                                   ///< Real raw noise, which will be merge into the driving force
        fftw_complex*                    fftw_In;                               ///< FFTW input pointer
        fftw_complex*                    fftw_Out;                              ///< FFTW out pointer
        fftw_plan                        FFTBackward;                           ///< Plan for the execution of FFTW
        std::complex<double>*            RandomFourier;                         ///< Pointer to storage of Fourier-space random number distribution
        std::complex<double>*            RandomReal;                            ///< Pointer to storage of real space random number distribution
        std::default_random_engine       generator;                             ///< Random number generator
        std::normal_distribution<double> distribution;                          ///< Random number distribution
    private:
};
}
#endif
