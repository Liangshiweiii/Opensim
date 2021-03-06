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

#include "Tools/Includes.h"
#include "DrivingForce.h"
#include "Info.h"
#include "Noise.h"
#include "Settings.h"
#include "Temperatures.h"
#include "VTK.h"
#include "fftw3.h"

namespace opensim
{
using namespace std;

Noise::Noise(Settings& locSettings)
{
    Initialize(locSettings);
    ReadInput(DefaultInputFileName);
}

Noise::Noise(Settings& locSettings, std::string InputFileName)
{
    Initialize(locSettings, InputFileName);
}

void Noise::Initialize(Settings& locSettings)
{
    thisclassname        = "Noise";
    //DefaultInputFileName = ProjectInputDir + "Noise.opi";

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    dx = locSettings.dx;
    dy = locSettings.dx;
    dz = locSettings.dx;

    dt = locSettings.dt;

    Nphases   = locSettings.Nphases;

    // Allocate three-dimensional Fourier space noise
    RandomFourier = new complex<double>[Nx*Ny*Nz]();
    RandomReal    = new complex<double>[Nx*Ny*Nz]();

    // Allocate raw Noise Storage
    Raw.Allocate(Nx, Ny, Nz, 0);

    // Assign storages to fftw pointers
    fftw_In  = reinterpret_cast<fftw_complex*>(RandomFourier);
    fftw_Out = reinterpret_cast<fftw_complex*>(RandomReal);

    // Create FFT-Plan, which determines how the FFT will be executed
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
    FFTBackward = fftw_plan_dft_3d(Nx,Ny,Nz,fftw_In,fftw_Out,FFTW_BACKWARD,FFTW_ESTIMATE);

    // Initialize random number distribution
    distribution = normal_distribution<double>(0.0,0.5);

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void Noise::Initialize(Settings& locSettings, std::string InputFileName)
{
    Initialize(locSettings);
    ReadInput(InputFileName);
}

void Noise::ReadInput(std::string InputFileName)
{
    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened",
                "ReadInput");
        exit(1);
    };

    Info::WriteBlankLine();
    Info::WriteLineInsert("Noise");
    Info::WriteStandard("Source", InputFileName);

    // Read Parameters
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
    
    gamma     = double(UserInterface::ReadParameterD(inp, moduleLocation, string("dgamma")));
    TimeSteps = int(   UserInterface::ReadParameterD(inp, moduleLocation, string("iTime")));
    KCutOff   = int(   UserInterface::ReadParameterD(inp, moduleLocation, string("iCutOff")));

    // Check if the input TimeSteps makes sense
    if (TimeSteps < 1)
    {
        std::string message = "A input value TimeSteps < 1 does not make sense! "
            "Set TimeSteps to 1!";
        Info::WriteWarning( message, thisclassname, "ReadInput");
        TimeSteps = 0;
    }

    // Check if the input KCutOff makes sense
    if (KCutOff < 0)
    {
        std::string message = "A input value KCutOff < 0 does not make sense! "
            "Set KCutOff to 0!";
        Info::WriteWarning( message, thisclassname, "ReadInput");
        KCutOff = 0;
    }
    else if (KCutOff > min(Nx,min(Ny,Nz)))
    {
        std::string message = "A input value for KCutOff bigger than the grid"
            " size does not make sense! Set KCutOff to 0";
        Info::WriteWarning( message, thisclassname, "ReadInput");
        KCutOff = 0;
    }
    Info::WriteLine();
    inp.close();
}

void Noise::Generate(void)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        if((i < Nx-KCutOff) and (j < Ny-KCutOff) and (k < Nz-KCutOff))
        {
            // Generate Gaussian white noise coefficients
            complex<double> cijk{distribution(generator),distribution(generator)};
            //cijk /= abs(cijk);

            // Store computed coefficient
            RandomFourier[k + Nz*(j + Ny*i)] = cijk;
        }
        else
        {
            RandomFourier[k + Nz*(j + Ny*i)] = 0.0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    // Do Fourier transform
    fftw_execute(FFTBackward);

    // Calculate global maximum of RandomReal
    double Max = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,reduction(max:Max))
    {
        Max = max(Max, abs(real(RandomReal[k + Nz*(j + Ny*i)])));
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        RandomReal[k + Nz*(j + Ny*i)] /= Max;
        RandomReal[k + Nz*(j + Ny*i)] -= Raw(i,j,k);
        RandomReal[k + Nz*(j + Ny*i)] /= TimeSteps;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Noise::AddToDrivingForce(int tStep, Temperature& Temp, DrivingForce& DF)
{
    if (tStep == 0)
    {
        // Generate initial noise distribution
        Generate();

        // Save initial noise distribution and add it to the driving force

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
        {
            // Set initial noise
            Raw(i,j,k) = TimeSteps * real(RandomReal[k + Nz*(j + Ny*i)]);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        // Generate new noise distribution in order to able to interpolate
        // between the old distribution and the new one.
        Generate();
    }
    else
    {
        // First check if a new noise Distribution has to be generated
        if (!(tStep%TimeSteps)) Generate();

        // Interpolate between old and new noise distribution by adding RandomReal
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
        {
            // Compute noise for this time step
            Raw(i,j,k) += real(RandomReal[k + Nz*(j + Ny*i)]);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }

    // Calculate the global maximum of Raw
    double Max = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,reduction(max:Max))
    {    
        Max = max(Max, abs(Raw(i,j,k)));
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    
    // Interpolate between old and new noise distribution by adding RandomReal
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        Raw(i,j,k) /= Max;

        // Add noise to driving force
        for (auto it = DF.Raw(i,j,k).cbegin(); it < DF.Raw(i,j,k).cend(); ++it)
        {
            double Ampl = sqrt( 2 * PhysicalConstants::kBoltzmann * Temp(i,j,k)/gamma);
            DF.Raw(i,j,k).add_asym(it->indexA, it->indexB, Ampl * Raw(i,j,k));
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void Noise::WriteVTK(int tStep)
{
    stringstream buffer;
    std::vector<int> DataTypes {PDScalars};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        VTK::WriteScalar(buffer, Raw, "Noise");
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "Noise", tStep);
}
}// end of name space
