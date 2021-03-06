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
#include "Tools/UserInterface.h"
#include "Tools/Node.h"
#include "Tools/NodeV.h"
#include "InterfaceEnergy.h"
#include "Settings.h"
#include "InterfaceDiffusion.h"
#include "PhaseField.h"
#include "VTK.h"
#include <exception>

namespace opensim
{
using namespace std;

InterfaceDiffusion::InterfaceDiffusion(const Settings& locSettings)
{
    Initialize(locSettings);
    ReadInput();
}

void InterfaceDiffusion::Initialize(const Settings& locSettings)
{
    thisclassname        = "InterfaceDiffusion";
    //DefaultInputFileName = ProjectInputDir + "InterfaceDiffusionInput.opi";

    Nphases = locSettings.Nphases;
    Nx      = locSettings.Nx;
    Ny      = locSettings.Ny;
    Nz      = locSettings.Nz;
    dx      = locSettings.dx;
    dy      = locSettings.dx;
    dz      = locSettings.dx;

    Coefficients.Allocate(Nphases, Nphases);

    DiffusionPotential.Allocate(Nx, Ny, Nz, 2);

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void InterfaceDiffusion::ReadInput(const std::string InputFileName)
{
    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened",
                "ReadInput");
        exit(1);
    };

    Info::WriteBlankLine();
    Info::WriteLineInsert("Interface diffusion");
    Info::WriteStandard("Source", InputFileName.c_str());

    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    DoubleObstacleSmoothnessRange =
        UserInterface::ReadParameterD(inp, moduleLocation, string("dSR"));

    for(unsigned int alpha = 0;     alpha < Nphases; alpha++)
    for(unsigned int beta  = alpha; beta  < Nphases;  beta++)
    {
        double value = 0.0;
        stringstream converter;
        converter << alpha << "_" << beta;
        string counter = converter.str();
        value = UserInterface::ReadParameterD(inp, moduleLocation, string("dIDC_") + counter,
                false, 0);
        Coefficients(alpha, beta) = value;
        Coefficients(beta, alpha) = value;
    }
    inp.close();
    Info::WriteBlankLine();
}

void InterfaceDiffusion::CalculatePhaseFieldIncrements(PhaseField& Phase,
        const InterfaceEnergy& Sigma)
{
    const unsigned int PhaseFieldBCells = Phase.Fields.Bcells();
    if (PhaseFieldBCells <= 2)
    {
        Info::WriteExit("Not enough PhaseField boundary cells (minimum 3 cells!)",
                thisclassname, "Calculate");
        exit(1);
    }

    const unsigned int InterfaceEnergyBCells = Sigma.Bcells();
    if (InterfaceEnergyBCells <= 1)
    {
        Info::WriteExit("Not enough InterfaceEnergy boundary cells (minimum 2 cells!)",
                thisclassname, "Calculate");
        exit(1);
    }

    // Override flag functionality
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields, 3,)
        Phase.Fields(i,j,k).flag = 2;
    OMP_PARALLEL_STORAGE_LOOP_END

    CalculateDiffusionPotential          (Phase, Sigma);
    CalculateDiffusionPotentialLaplacian (Phase);
}

void InterfaceDiffusion::CalculateDiffusionPotentialLaplacian(
        PhaseField& Phase)
{
    const double dx_2 = 1.0/(dx*dx);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DiffusionPotential, 0,)
    for (int ii = -1; ii <= +1; ++ii)
    for (int jj = -1; jj <= +1; ++jj)
    for (int kk = -1; kk <= +1; ++kk)
    {
        for (auto it = DiffusionPotential(i+ii,j+jj,k+kk).cbegin();
                  it < DiffusionPotential(i+ii,j+jj,k+kk).cend(); ++it)
        {
            const int pIndexA = Phase.FieldsStatistics[it->indexA].Phase;
            const int pIndexB = Phase.FieldsStatistics[it->indexB].Phase;

            double value = - Coefficients(pIndexA,pIndexB)
                * LaplacianStencil27[ii+1][jj+1][kk+1]*(it->value)*dx_2;
            Phase.FieldsDot(i,j,k).add(it->indexA, it->indexB, value);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DiffusionPotential, 1,)
        DiffusionPotential(i,j,k).clear();
    OMP_PARALLEL_STORAGE_LOOP_END
}

double InterfaceDiffusion::PotentialDerivative(const double alpha,
        const double beta) const
{
    double value = 0.0;
    // This method calulates the
    if (DoubleObstacleSmoothnessRange == 0.)
    {
        // Normal double-obstacle potential
        double Pot = alpha * beta;
        if (Pot < 0.)
            value = - beta;
        else if (Pot > 0.)
            value = beta;
    }
    else if (DoubleObstacleSmoothnessRange == 1.)
    {
        // Normal double-well potential
        value =  2 * alpha * pow(beta,2);
    }
    else
    {
        // Bent-cable model for the potential
        // Calculate bent-cable function of phi-beta
        double BentCableBeta = 0.0;
        if (abs(beta) < DoubleObstacleSmoothnessRange)
        {
            BentCableBeta -= 1.0 * pow(beta,4)
                / (16.0 * pow(DoubleObstacleSmoothnessRange,3));
            BentCableBeta += 3.0 * pow(beta,2)
                / ( 8.0 * pow(DoubleObstacleSmoothnessRange,1));
            BentCableBeta += 1.0 * pow(beta,1) /   2.0;
            BentCableBeta += 3.0 * DoubleObstacleSmoothnessRange / 16.0;
        }
        else if (beta >= DoubleObstacleSmoothnessRange) BentCableBeta = beta;

        // Calculate smoothed absolute value of phi-alpha
        double SAbsBeta  = - beta  + 2.0 * BentCableBeta;

        // Calculate bent-cable function of phi-alpha
        double dBentCableAlpha_dAlpha = 0.0;
        if (abs(alpha) < DoubleObstacleSmoothnessRange)
        {
            dBentCableAlpha_dAlpha -= 1.0 * pow(alpha,3)
                / (4.0 * pow(DoubleObstacleSmoothnessRange,3));
            dBentCableAlpha_dAlpha += 6.0 * pow(alpha,1)
                / (8.0 * pow(DoubleObstacleSmoothnessRange,1));
            dBentCableAlpha_dAlpha += 0.5;
        }
        else if (alpha >= DoubleObstacleSmoothnessRange)
            dBentCableAlpha_dAlpha = 1.0;

        // Calculate smoothed absolute value of phi-alpha
        double dSAbsAlpha_dAlpha = - 1.0 + 2.0 * dBentCableAlpha_dAlpha;

        return dSAbsAlpha_dAlpha * SAbsBeta;
    }
    return value;
}

void InterfaceDiffusion::CalculateDiffusionPotential(PhaseField& Phase,
        const InterfaceEnergy& Sigma)
{
    const double Prefactor2 = Pi*Pi/(Phase.Eta*Phase.Eta);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields, 2,)
    {
        double norm_1 = 1.0/double(Phase.Fields(i,j,k).size());

        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha < Phase.Fields(i,j,k).cend() - 1; ++alpha)
        for(auto  beta = alpha + 1;
                  beta < Phase.Fields(i,j,k).cend(); ++beta)
        {
            double dDiffusionPotential_dt =
                Sigma(i, j, k, alpha->index, beta->index) *
                ((alpha->laplacian + Prefactor2 *
                  PotentialDerivative(beta->value,  alpha->value)) -
                 (beta->laplacian + Prefactor2 *
                  PotentialDerivative(alpha->value, beta->value)));

            if(Phase.Fields(i,j,k).size() > 2)
            for(auto gamma = Phase.Fields(i,j,k).cbegin();
                     gamma < Phase.Fields(i,j,k).cend(); ++gamma)
            if((gamma != alpha) && (gamma != beta))
            {
                dDiffusionPotential_dt +=
                    Sigma(i, j, k, beta->index,  gamma->index) *
                    (gamma->laplacian + Prefactor2 *
                     PotentialDerivative(beta->value,  gamma->value)) -
                    Sigma(i, j, k, alpha->index, gamma->index) *
                    (gamma->laplacian + Prefactor2 *
                     PotentialDerivative(alpha->value, gamma->value));
            }

            dDiffusionPotential_dt *= norm_1;

            DiffusionPotential(i,j,k).add_asym(alpha->index, beta->index,
                    dDiffusionPotential_dt);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

InterfaceDiffusionAnisotropic::InterfaceDiffusionAnisotropic(const Settings& locSettings)
{
    Initialize(locSettings);
    InterfaceDiffusion::ReadInput(DefaultInputFileName);

}

void InterfaceDiffusionAnisotropic::Initialize(const Settings& locSettings)
{
    InterfaceDiffusion::Initialize(locSettings);
    DiffusionFlux.Allocate     (Nx, Ny, Nz, 1);
}

void InterfaceDiffusionAnisotropic::ReadInput(const std::string InputFileName)
{
    InterfaceDiffusion::ReadInput(InputFileName);
    thisclassname = "InterfaceDiffusionAnisotropic";
    if (DoubleObstacleSmoothnessRange < 0.04)
    {
        Info::WriteWarning("DoubleObstacleSmoothnessRange might be too small to the algorithm to work!",
                thisclassname, "ReadInput");
    }
}

void InterfaceDiffusionAnisotropic::CalculatePhaseFieldIncrements(PhaseField& Phase,
        const InterfaceEnergy& Sigma)
{
    const int PhaseFieldBCells = Phase.Fields.Bcells();
    if (PhaseFieldBCells <= 2)
    {
        Info::WriteExit("Not enough PhaseField boundary cell (min 3 cells!)",
                thisclassname, "CalculateAnsiotropic");
        exit(1);
    }

    const int InterfaceEnergyBCells = Sigma.Bcells();
    if (InterfaceEnergyBCells <= 1)
    {
        Info::WriteExit("Not enough InterfaceEnergy boundary cell (min 2 cells!)",
                thisclassname, "CalculateAnsiotropic");
        exit(1);
    }

    // Override flag functionality
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields, 3,)
        Phase.Fields(i,j,k).flag = 2;
    OMP_PARALLEL_STORAGE_LOOP_END

    // This method solves the anisotropic diffusion-equation
    // PhiDot = Nabla ( Diffusion-Tensor * Nabla( Diffusion-potential))
    CalculateDiffusionPotential          (Phase, Sigma);
    CalculateDiffusionPotentialGradients ();
    CalculateDiffusionFlux               (Phase);
    CalculateDiffusionFluxDivergence     (Phase);
}

void InterfaceDiffusionAnisotropic::CalculateDiffusionFlux(PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DiffusionFlux, 1,)
        NodeV locNormals       = Phase.Normals(i,j,k);
        NodeV locGradients     = Phase.Gradients(i,j,k);
        NodeV locDiffusionFlux = DiffusionFlux(i,j,k);

        for (auto it = locDiffusionFlux.cbegin();
                  it < locDiffusionFlux.cend(); ++ it)
        {
            dMatrix3x3 locDT;
            locDT.set_to_unity();

            // Calculation of the anisotropy:
            // The ratio of tangential to normal Interface Diffusion will be
            // calculated
            //TODO Note: The bent-cable potential is not used !!
            double Ratio     = 0.0; // Ration of normal to tangential diffusion
            double Potential = 0.0; // Stores value of double obstacle potential
            Potential  = abs(Phase.Fields(i,j,k).get(it->indexA));
            Potential *= abs(Phase.Fields(i,j,k).get(it->indexB));
            if (Potential > 1.0e-200)
            {
               Ratio  = abs(locGradients.get(it->indexA)
                       * locGradients.get(it->indexB));
               Ratio /= Potential;
               Ratio *= (Phase.Eta * Phase.Eta)/(Pi * Pi);
            }

            // Limit Ratio to meaningful values
            if (Ratio > 1.0) Ratio = 1.0; // Extinguishes normal diffusion completely
            if (Ratio < 0.0) Ratio = 0.0; // Results in the scalar model

            // Calculate local interface diffusion mobility tensor
            dVector3 locNormal = locNormals.get(it->indexA, it->indexB);
            locDT(0,0) -= Ratio * locNormal.getX()*locNormal.getX();
            locDT(0,1) -= Ratio * locNormal.getX()*locNormal.getY();
            locDT(0,2) -= Ratio * locNormal.getX()*locNormal.getZ();
            locDT(1,0) -= Ratio * locNormal.getY()*locNormal.getX();
            locDT(1,1) -= Ratio * locNormal.getY()*locNormal.getY();
            locDT(1,2) -= Ratio * locNormal.getY()*locNormal.getZ();
            locDT(2,0) -= Ratio * locNormal.getZ()*locNormal.getX();
            locDT(2,1) -= Ratio * locNormal.getZ()*locNormal.getY();
            locDT(2,2) -= Ratio * locNormal.getZ()*locNormal.getZ();

            // Apply diffuse interface interpolation
            int pIndexA = Phase.FieldsStatistics[it->indexA].Phase;
            int pIndexB = Phase.FieldsStatistics[it->indexB].Phase;
            locDT *= - Coefficients(pIndexA,pIndexB);

            // Apply Projection Operator
            dVector3 value = locDT * locDiffusionFlux.get(it->indexA,it->indexB);
            DiffusionFlux(i,j,k).set(it->indexA, it->indexB, value);
        }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InterfaceDiffusionAnisotropic::CalculateDiffusionFluxDivergence(PhaseField& Phase)
{
    const double DWeights[3] = {-0.5/dx, 0.0, 0.5/dx};

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DiffusionFlux, 0,)
    for (int ii = -1; ii <= +1; ii+=2)
    {
        double DWeight = DWeights[ii+1];
        for (auto it = DiffusionFlux(i+ii,j,k).cbegin();
                  it < DiffusionFlux(i+ii,j,k).cend(); ++it)
        {
            double valueX = DWeight*(it->X);
            Phase.FieldsDot(i,j,k).add_asym(it->indexA, it->indexB, valueX);
        }
        for (auto it = DiffusionFlux(i,j+ii,k).cbegin();
                  it < DiffusionFlux(i,j+ii,k).cend(); ++it)
        {
            double valueY = DWeight*(it->Y);
            Phase.FieldsDot(i,j,k).add_asym(it->indexA, it->indexB, valueY);
        }
        for (auto it = DiffusionFlux(i,j,k+ii).cbegin();
                  it < DiffusionFlux(i,j,k+ii).cend(); ++it)
        {
            double valueZ = DWeight*(it->Z);
            Phase.FieldsDot(i,j,k).add_asym(it->indexA, it->indexB, valueZ);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    const int BCells = DiffusionFlux.Bcells();
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DiffusionFlux, BCells,)
        DiffusionFlux(i,j,k).clear();
    OMP_PARALLEL_STORAGE_LOOP_END
}

void InterfaceDiffusionAnisotropic::CalculateDiffusionPotentialGradients()
{
    const double GWeights[3] = {-0.5/dx, 0.0, 0.5/dx};

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DiffusionPotential, 1,)
        DiffusionFlux(i,j,k).clear();
        for (int ii = -1; ii <= +1; ii+=2)
        {
            double GWeight = GWeights[ii+1];
            for (auto it = DiffusionPotential(i+ii,j,k).cbegin();
                      it < DiffusionPotential(i+ii,j,k).cend(); ++it)
            {
                double value = GWeight*(it->value);
                DiffusionFlux(i,j,k).add_X(it->indexA, it->indexB, value);
            }
            for (auto it = DiffusionPotential(i,j+ii,k).cbegin();
                      it < DiffusionPotential(i,j+ii,k).cend(); ++it)
            {
                double value = GWeight*(it->value);
                DiffusionFlux(i,j,k).add_Y(it->indexA, it->indexB, value);
            }
            for (auto it = DiffusionPotential(i,j,k+ii).cbegin();
                      it < DiffusionPotential(i,j,k+ii).cend(); ++it)
            {
                double value = GWeight*(it->value);
                DiffusionFlux(i,j,k).add_Z(it->indexA, it->indexB, value);
            }
        }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,DiffusionPotential, 2,)
        DiffusionPotential(i,j,k).clear();
    OMP_PARALLEL_STORAGE_LOOP_END
}

}
