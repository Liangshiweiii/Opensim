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

#include "HeatDiffusion.h"
#include "Temperatures.h"
#include "Info.h"
#include "Settings.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"

namespace opensim
{
using namespace std;

HeatDiffusion::HeatDiffusion(const Settings& locSettings)
{
    Initialize(locSettings);
    ReadInput();
}

HeatDiffusion::HeatDiffusion(const Settings& locSettings, const std::string FileName)
{
    Initialize(locSettings);
    ReadInput(FileName);
}

void HeatDiffusion::Initialize(const Settings& locSettings)
{
    thisclassname = "HeatDiffusion";
    //DefaultInputFileName = ProjectInputDir + "HeatDiffusionInput.opi";

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;
    Nphases = locSettings.Nphases;
    dx = locSettings.dx;
    dx2Inv = 1.0/(dx*dx);

    PhaseThermalDiffusivity.Allocate(Nphases);
    PhaseHeatCapacity.Allocate(Nphases);
    //PhaseDensity.Allocate(Nphases); //NOT USED

    EffectiveThermalDiffusivity.Allocate(Nx, Ny, Nz, 0);
    EffectiveHeatCapacity.Allocate(Nx, Ny, Nz, 0);
    //EffectiveDensity.Allocate(Nx, Ny, Nz, 0); //NOT USED

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void HeatDiffusion::Initialize(const Settings& locSettings, const std::string FileName)
{
    Initialize(locSettings);
    if (FileName.empty() or FileName == "default") ReadInput();
    else ReadInput(FileName);
}

void HeatDiffusion::ReadInput(const std::string FileName)
{
    fstream inp(FileName.c_str(), ios::in);

    if (!inp)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };

    Info::WriteBlankLine();
    Info::WriteLineInsert("Heat diffusion properties");
    Info::WriteStandard("Source", FileName.c_str());
    
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
    
    for(int n = 0; n < Nphases; n++)
    {
        stringstream converter;
        converter << n;
        string counter = converter.str();
        Info::WriteBlankLine();
        Info::WriteStandard("Heat diffusion properties for phase", std::to_string(n));
        PhaseThermalDiffusivity[n] = UserInterface::ReadParameterD(inp, moduleLocation, string("ThermalDiffusivity_") + counter);
        PhaseHeatCapacity[n] = UserInterface::ReadParameterD(inp, moduleLocation, string("HeatCapacity_") + counter);
        //PhaseDensity[n] = UserInterface::ReadParameterD(inp, string("Density_") + counter, true, 0);
    }
    Info::WriteLine();
    inp.close();

    MaxThermalDiffusivity = 0.0;
    for(int n = 0; n < Nphases; n++)
    {
        if (MaxThermalDiffusivity < PhaseThermalDiffusivity[n])
        {
            MaxThermalDiffusivity = PhaseThermalDiffusivity[n];
        }
    }
}

void HeatDiffusion::Remesh(const int newNx, const int newNy, const int newNz)
{
    EffectiveThermalDiffusivity.Reallocate(newNx, newNy, newNz);
    EffectiveHeatCapacity.Reallocate(newNx, newNy, newNz);
    //EffectiveDensity.Reallocate(newNx,newNy, newNz);

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    Info::WriteStandard(thisclassname, "Remeshed");
}

void HeatDiffusion::SetEffectiveProperties(const PhaseField& Phase)
{
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    {
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i,j,k).cbegin(); alpha < Phase.Fields(i,j,k).cend(); ++alpha)
            {
                int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                EffectiveHeatCapacity(i,j,k) += PhaseHeatCapacity[pIndexA]*alpha -> value;
                EffectiveThermalDiffusivity(i,j,k) += PhaseThermalDiffusivity[pIndexA]*alpha -> value;
                //EffectiveDensity(i,j,k) += PhaseDensity[pIndexA]*alpha -> value;
            }
        }
        else
        {
            int pIndexA = Phase.FieldsStatistics[Phase.Fields(i,j,k).front().index].Phase;
            EffectiveHeatCapacity(i,j,k) = PhaseHeatCapacity[pIndexA];
            EffectiveThermalDiffusivity(i,j,k) = PhaseThermalDiffusivity[pIndexA];
            //EffectiveDensity(i,j,k) = PhaseDensity[pIndexA];
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void HeatDiffusion::CalculateLaplacian(Temperature& Temperature) const
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Temperature.TxDx,0,)
    {
        Temperature.TxDx(i,j,k)[0] = dx2Inv*
                        (Temperature.Tx(i+1,j,k) - 2.0*Temperature.Tx(i,j,k) + Temperature.Tx(i-1,j,k));

        Temperature.TxDx(i,j,k)[1] = dx2Inv*
                        (Temperature.Tx(i,j+1,k) - 2.0*Temperature.Tx(i,j,k) + Temperature.Tx(i,j-1,k));

        Temperature.TxDx(i,j,k)[2] = dx2Inv*
                        (Temperature.Tx(i,j,k+1) - 2.0*Temperature.Tx(i,j,k) + Temperature.Tx(i,j,k-1));
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void HeatDiffusion::UpdateTemperature(Temperature& Temperature, const double dt) const
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Temperature.Tx,0,)
    {
        Temperature.Tx(i,j,k) += dt*Temperature.qdot(i,j,k)/EffectiveHeatCapacity(i,j,k) +
            dt*EffectiveThermalDiffusivity(i,j,k)*
            (Temperature.TxDx(i,j,k)[0] + Temperature.TxDx(i,j,k)[1] + Temperature.TxDx(i,j,k)[2]);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    //SetBoundaryConditions();
}

void HeatDiffusion::ClearEffectiveProperties()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EffectiveHeatCapacity,0,)
    {
        EffectiveThermalDiffusivity(i,j,k) = 0.0;
        EffectiveHeatCapacity(i,j,k) = 0.0;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void HeatDiffusion::Solve(const PhaseField& Phase, Temperature& Temperature,
        const BoundaryConditions& BC, const double dt)
{
    // Check stability

    double eta = 0.25*(dx*dx)/MaxThermalDiffusivity;

    if(dt > eta)
    {
        std::string message = "Stability criterion.\n dt (" +
                to_string(dt) + ") > eta (" + to_string(eta) + ")\n";
        Info::WriteExit(thisclassname, message, "Solve()");
        exit(112);
    }

    // Solve

    Temperature.SetBoundaryConditions(BC);
    SetEffectiveProperties(Phase);
    CalculateLaplacian(Temperature);
    UpdateTemperature(Temperature, dt);
    ClearEffectiveProperties();
}

}// namespace openphase
