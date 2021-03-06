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
#include "Info.h"
#include "PhaseField.h"
#include "Orientations.h"
#include "Settings.h"
#include "BoundaryConditions.h"
#include "Temperatures.h"
#include "Tools/UserInterface.h"
#include "VTK.h"
#include "Velocities.h"

namespace opensim
{
using namespace std;

PlasticProperties::PlasticProperties(Settings& locSettings)
{
    Initialize(locSettings);
}

PlasticProperties::~PlasticProperties()
{
}

void PlasticProperties::Initialize(Settings& locSettings)
{
    thisclassname = "PlasticProperties";
    //DefaultInputFileName = ProjectInputDir + "PlasticFlowInput.opi";

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    PlasticStrain.Allocate(Nx, Ny, Nz, 1);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,PlasticStrain,PlasticStrain.Bcells(),)
    {
        PlasticStrain(i,j,k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Info::WriteStandard(thisclassname, "Initialized");
}

void PlasticProperties::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(PlasticStrain);
    BC.SetY(PlasticStrain);
    BC.SetZ(PlasticStrain);
}

void PlasticProperties::Remesh(int newNx, int newNy, int newNz,
                                                         BoundaryConditions& BC)
{
    SetBoundaryConditions(BC);

    PlasticStrain.Remesh(newNx, newNy, newNz);

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    SetBoundaryConditions(BC);
}

void PlasticProperties::Write(int tStep, bool legacy_format)
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"PlasticStrain_", tStep, ".dat");

    ofstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be created", thisclassname, "Write()");
        exit(1);
    };

    if(not legacy_format)
    {
        int tmp = Nx;
        out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
        tmp = Ny;
        out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
        tmp = Nz;
        out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
    }
    STORAGE_LOOP_BEGIN(i,j,k,PlasticStrain,0)
    {
        for (int n = 0; n < 6; n++)
        {
            double temp = PlasticStrain(i,j,k)[n];
            out.write(reinterpret_cast<char*>(&temp), sizeof(double));
        }
    }
    STORAGE_LOOP_END
    out.close();
}

void PlasticProperties::Read(const BoundaryConditions& BC, int tStep, bool legacy_format)
{
    string FileName =
        UserInterface::MakeFileName(RawDataDir,"PlasticStrain_", tStep, ".dat");

    Read(FileName, legacy_format);
    SetBoundaryConditions(BC);
}

void PlasticProperties::Read(string FileName, bool legacy_format)
{
    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit(FileName + " could not be opened",
                thisclassname, "Read()");
        exit(1);
    };

    if(not legacy_format)
    {
        int locNx = Nx;
        int locNy = Ny;
        int locNz = Nz;
        inp.read(reinterpret_cast<char*>(&locNx), sizeof(int));
        inp.read(reinterpret_cast<char*>(&locNy), sizeof(int));
        inp.read(reinterpret_cast<char*>(&locNz), sizeof(int));
        if(locNx != Nx or locNy != Ny or locNz != Nz)
        {
            stringstream message;
            message << "Inconsistent system dimensions!\n"
                    << "Input data dimensions: (" << locNx
                    << ", " << locNy << ", " << locNz << ") grid points.\n"
                    << "Required data dimensions: (" << Nx
                    << ", " << Ny << ", " << Nz << ") grid points.\n";
            Info::WriteExit(message.str(), thisclassname, "Read()");
            exit(1);
        }
    }
    STORAGE_LOOP_BEGIN(i,j,k,PlasticStrain,0)
    {
        double val = 0.0;
        PlasticStrain(i,j,k).set_to_zero();
        for (int n = 0; n < 6; n++)
        {
            inp.read(reinterpret_cast<char*>(&val), sizeof(double));            // Fields(i,j,k)->value
            PlasticStrain(i,j,k)[n] = val;
        }
    }
    STORAGE_LOOP_END
    inp.close();

    Info::WriteStandard(thisclassname, "Binary input loaded");
}

void PlasticProperties::Advect(Velocities& Vel, BoundaryConditions& BC,
                    double dt, int scheme)
{
    if(StrainsDot.IsNotAllocated())
    {
        StrainsDot.Allocate(Nx, Ny, Nz, 1);
    }
    if(not StrainsDot.IsSize(Nx, Ny, Nz))
    {
        StrainsDot.Reallocate(Nx, Ny, Nz);
    }
    SetBoundaryConditions(BC);

    switch(scheme)
    {
        case Upwind:
        {
            const double dx = Vel.dx;
            const double dx2 = 0.5/dx;
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,StrainsDot,0,)
            for (int n = 0; n < 6; n++)
            {
                StrainsDot(i,j,k)[n] = dx2 *
                     ((fabs(Vel.Average(i-1,j,k)[0]) + Vel.Average(i-1,j,k)[0])*PlasticStrain(i-1,j,k)[n] +
                      (fabs(Vel.Average(i,j-1,k)[1]) + Vel.Average(i,j-1,k)[1])*PlasticStrain(i,j-1,k)[n] +
                      (fabs(Vel.Average(i,j,k-1)[2]) + Vel.Average(i,j,k-1)[2])*PlasticStrain(i,j,k-1)[n] +
                      (fabs(Vel.Average(i+1,j,k)[0]) - Vel.Average(i+1,j,k)[0])*PlasticStrain(i+1,j,k)[n] +
                      (fabs(Vel.Average(i,j+1,k)[1]) - Vel.Average(i,j+1,k)[1])*PlasticStrain(i,j+1,k)[n] +
                      (fabs(Vel.Average(i,j,k+1)[2]) - Vel.Average(i,j,k+1)[2])*PlasticStrain(i,j,k+1)[n]) -
                      (fabs(Vel.Average(i,j,k)[0]) +
                       fabs(Vel.Average(i,j,k)[1]) +
                       fabs(Vel.Average(i,j,k)[2])) * PlasticStrain(i, j, k)[n]/dx;
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            break;
        }
        case LaxWendroff:
        {
            Info::WriteExit("LaxWendroff advection scheme not supported for PlasticStrain",
                            thisclassname, "AdvectPlasticStrain(PFP, Vel, BC, dt)");
            exit(13);
            break;
        }
        default:
        {
            Info::WriteExit("Wrong/None advection scheme given in the input file",
                             thisclassname, "AdvectPlasticStrain(PFP, Vel, BC, dt)");
            exit(13);
        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,StrainsDot,0,)
    {
        for (int n = 0; n < 6; ++n)
        {
            PlasticStrain(i, j, k)[n] += StrainsDot(i, j, k)[n]*dt;
            StrainsDot(i, j, k)[n] = 0.0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    SetBoundaryConditions(BC);
}

void PlasticProperties::WritePlasticStrainVTK(int tStep)
{
    stringstream outbufer;

    outbufer << "# vtk DataFile Version 3.0\n";
    outbufer << "PlasticStrain\n";
    outbufer << "ASCII\n";
    outbufer << "DATASET STRUCTURED_GRID\n";
    outbufer << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    outbufer << "POINTS " <<  Nx*Ny*Nz << " double\n";

    for(double k = 0; k < Nz; ++k)
    for(double j = 0; j < Ny; ++j)
    for(double i = 0; i < Nx; ++i)
    {
        outbufer << i << " " << j << " " << k << "\n";
    }
    outbufer << "\n";
    outbufer << "POINT_DATA " << Nx*Ny*Nz << "\n";

    for (int n = 0; n < 6; ++n)
    {
        outbufer << "SCALARS PlasticStrain_" << n << " double 1\n";
        outbufer << "LOOKUP_TABLE default\n";

        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            outbufer << PlasticStrain(i,j,k)[n] << "\n";
        }
    }

    outbufer << "SCALARS EquivalentPlasticStrain" << " double 1\n";
    outbufer << "LOOKUP_TABLE default\n";

    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        outbufer << PlasticStrain(i,j,k).norm() << "\n";
    }

    string FileName = UserInterface::MakeFileName(VTKDir, "PlasticStrain_",
                                                                 tStep, ".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbufer.rdbuf();
    vtk_file.close();
}

}// namespace openphase
