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

#ifndef PLASTICFLOWCP_H
#define PLASTICFLOWCP_H

#include "Tools/Includes.h"
#include "Tools/NodeVn.h"
#include "Mechanics/PlasticFlow/PlasticFlow.h"

namespace opensim
{
class Settings;
class PhaseField;
class ElasticProperties;
class Composition;
class Temperature;
class Settings;
class Orientations;
class BoundaryConditions;
class DrivingForce;
class Velocities;
class ChemicalProperties;
class Composition;

class PlasticFlowCP : public PlasticFlow
{
public:
    PlasticFlowCP(){};
    PlasticFlowCP(Settings& locSettings);
    ~PlasticFlowCP(){};

    void Initialize(Settings& locSettings);
    using OPObject::ReadInput;
    void ReadInput(std::string InputFileName);

    void SetBoundaryConditions(const BoundaryConditions& BC);
    void Remesh(int newNx, int newNy, int newNz, BoundaryConditions& BC);
    void Advect(PhaseField& Phase, Velocities& Vel, BoundaryConditions& BC,
                        double dt, int scheme = Upwind);
    void SetInitialHardening(const PhaseField& Phase, BoundaryConditions& BC);
    void SetInitialHardening(const PhaseField& Phase, BoundaryConditions& BC, int targetPhaseIndex);
    void SetInitialHardening(const PhaseField& Phase, const Composition& Cx, BoundaryConditions& BC);
    dVector6 GetSchmidMatrixSym(int i, int j, int k, int phaseIndex, int slipSys, PhaseField& Phase, ElasticProperties& EP, Orientations& OR);
    dMatrix3x3 GetSchmidMatrix(int i, int j, int k, int phaseIndex, int slipSys, PhaseField& Phase, ElasticProperties& EP, Orientations& OR);

    double Solve(std::vector<OPObject*> Objects, bool verbose);

    void Write(int tStep, bool legacy_format = true);
    void Read(const BoundaryConditions& BC, int tStep, bool legacy_format = true);
    void Read(std::string FileName, bool legacy_format);

    void WriteCRSSVTK(PhaseField& Phase, int tStep);

    int Nx;
    int Ny;
    int Nz;
    int Nphases;
    int nSlipSystems;
    int PlasticityStartStep;
    double dt;
    double allowedShearRate;
    double TotalCRSSHardening;
    bool FirstStepHardening;

    /*
     * first index: thermodynamic phase
     * second index: slip system
     */
    bool PlasticitySwitch;
    bool HardeningSwitch;
    Tensor <bool, 1> Hardening;
    Tensor <bool, 1> PlasticitySwitchPhase;
    Tensor <double, 1> phaseExponent;
    Tensor <double, 1> phaseSlipRate;
    Tensor <double, 1> phaseHardExponent;
    Tensor <double, 1> BurgersLength;
    Tensor <int, 1> Nslip;
    Tensor<double, 2> PhaseCriticalResolvedShearStress;
    Tensor<double, 2> InitialHardeningModulus;
    Tensor<double, 2> LatentHardeningParameter;
    Tensor<double, 2> PhaseSaturationCRSS;
//    Tensor<dMatrix3x3, 2> PhaseSchmidMatrices;
    Storage3D<vStrain,0> MaxPlasticStrain;
    Storage<Tensor<double, 2>> PhaseGSNormal;                                   /// Storage: Matrix of glide plane normals for FCC (12 GS a 3 vector components)
    Storage<Tensor<double, 2>> PhaseGSDirection;                                /// Storage: Matrix of glide directions for FCC (12 GS a 3 vector components)
    Storage<Tensor<double, 2>> PhaseGSLine;                                     /// Storage: Matrix of dislocation lines for FCC (12 GS a 3 vector components)

    Storage3D<NodeVn<12>, 0> CRSS;
    Storage3D<NodeVn<12>, 0> CRSSdot;                                           /// storage for the increment of the CRSS due to advection

protected:
private:

};
}// namespace openphase

#endif

