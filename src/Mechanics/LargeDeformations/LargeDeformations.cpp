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

#include "Mechanics/Storages/PlasticProperties.h"
#include "Mechanics/LargeDeformations/LargeDeformations.h"
#include "Mechanics/LargeDeformations/LargeDeformationsMethods.h"
#include "Info.h"
#include "Velocities.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "Orientations.h"
#include "Settings.h"
#include "DrivingForce.h"
#include "PhaseField.h"
#include "Temperatures.h"
#include "Mechanics/ElasticitySolvers/ElasticitySolverSpectral.h"
#include "Compositions.h"
#include "BoundaryConditions.h"
#include "Tools/Quaternion.h"
#include "Tools/UserInterface.h"
#include "EquilibriumPartitionDiffusion.h"
#include "InterfaceEnergy.h"
#include "InterfaceMobility.h"
#include "Mechanics/Storages/DamageProperties.h"
#include "Mechanics/DamageModels/DamageModel.h"
#include "Mechanics/PlasticFlow/PlasticFlow.h"
#include "Mechanics/PlasticFlow/PlasticFlowCP.h"
#include "Mechanics/PlasticFlow/PlasticFlowMethods.h"
#include "DoubleObstacle.h"

namespace opensim
{
using namespace std;

LargeDeformations::LargeDeformations(Settings& locSettings,
                                                          ElasticProperties& EP)
{
    this->Initialize(locSettings, EP);
}

void LargeDeformations::Initialize(Settings& locSettings,
                                                          ElasticProperties& EP)
{
    thisclassname = "LargeDeformations";
    //DefaultInputFileName = ProjectInputDir + "LDInput.opi";
    dx = locSettings.dx;
	dt = locSettings.dt;

    counter = 0;
    maxIterations = 100;

	EP.InitializeLD();

    MaxAllowedStrainIncrement = 0.01;
    StrainAccuracy = 1.0e-5;                                                    // Strain accuracy
    MaxPressure = 100.0;                                                        // Maximal allowed pressure
    MaxSolverIterations = 200;

    IntfStabilisationNeeded = false;

//    PhaseDamageEnergy.resize(5000, 0.0);
//    PhasePlasticEnergy.resize(5000, 0.0);

    EffectiveEigenStrainsOLD.Allocate(EP.Nx, EP.Ny, EP.Nz, 1);
    RotationsOLD.Allocate(EP.Nx, EP.Ny, EP.Nz, 1);
    StiffnessOLD.Allocate(EP.Nx, EP.Ny, EP.Nz, 1);
    CRSSOLD.Allocate(EP.Nx, EP.Ny, EP.Nz, 1);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, RotationsOLD, RotationsOLD.Bcells(),)
    {
        RotationsOLD(i,j,k).set_to_unity();
        StiffnessOLD(i,j,k).set_to_zero();
        for(int ss = 0; ss < 12; ss++)
        {
            CRSSOLD(i,j,k).set_to_zero(ss);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    
    VelocityBC.Initialize(locSettings);
    VelocityBC.BC0X = 2;
    VelocityBC.BCNX = 2;
    VelocityBC.BC0Y = 2;
    VelocityBC.BCNY = 2;
    VelocityBC.BC0Z = 2;
    VelocityBC.BCNZ = 2;
//    if(EP.VelocityGradient.IsEmpty()) // Allocation moved back to EP.
//    {
//        EP.VelocityGradient.Allocate(EP.Nx, EP.Ny, EP.Nz, 1);
//    }

    initialized = true;
    Info::WriteLine();
    Info::WriteStandard(thisclassname, "Initialized");
}

void LargeDeformations::ReadInput(string InputFileName)
{
    Info::WriteLine();

    fstream inp(InputFileName.c_str(), ios::in);

    if (!inp or DefaultInputFileName.empty())
    {
        Info::WriteWarning("File " + InputFileName + " could not be opened\n" +
            "          Using standard parameters" ,thisclassname, "Initialize");
    }
    else
    {
		Info::WriteLineInsert("Large Deformations");
		Info::WriteStandard("Source", InputFileName);
		int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
        Verbose = UserInterface::ReadParameterB(inp, moduleLocation,
                string("verbose"), false, "No");
		MaxSolverIterations = UserInterface::ReadParameterI(inp, moduleLocation, string("MaxSolverIterations"),
			false, MaxSolverIterations);
        StrainAccuracy = UserInterface::ReadParameterD(inp, moduleLocation, string("SolverPrec1"), false, StrainAccuracy);
        MaxPressure = UserInterface::ReadParameterD(inp, moduleLocation, string("SolverPrec2"), false, MaxPressure);
        MaxAllowedStrainIncrement = UserInterface::ReadParameterD(inp, moduleLocation, string("MaxStrainInc"),
                                                  false, MaxAllowedStrainIncrement);
        IntfStabilisationNeeded = UserInterface::ReadParameterB(inp, moduleLocation,
                              string("IntfStabilisationNeeded"), false, "No");
        minIterations = UserInterface::ReadParameterI(inp, moduleLocation, string("MinLDIterations"),
                                          false, 4);
		maxIterations = UserInterface::ReadParameterI(inp, moduleLocation, string("MaxLDIterations"),
			false, 32);
    }
	Info::WriteLine();
}

LargeDeformations::~LargeDeformations(void)
{
}

void LargeDeformations::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(EffectiveEigenStrainsOLD);
    BC.SetY(EffectiveEigenStrainsOLD);
    BC.SetZ(EffectiveEigenStrainsOLD);

    BC.SetX(RotationsOLD);
    BC.SetY(RotationsOLD);
    BC.SetZ(RotationsOLD);

    BC.SetX(StiffnessOLD);
    BC.SetY(StiffnessOLD);
    BC.SetZ(StiffnessOLD);
}

void LargeDeformations::Remesh(int Nx, int Ny, int Nz,
                                                         BoundaryConditions& BC)
{
    EffectiveEigenStrainsOLD.Remesh(Nx, Ny, Nz);
    StiffnessOLD.Remesh(Nx, Ny, Nz);

    //Remeshing of RotationsOLD:
    Storage3D<Quaternion, 0>    QuaternionsOLD;
//    cout<<"new: Nx = "<<Nx<<"; Ny = "<<Ny<<"; Nz = "<<Nz<<endl;
//    cout<<"old: Nx = "<<RotationsOLD.sizeX()<<"; Ny = "<<RotationsOLD.sizeY()<<"; Nz = "<<RotationsOLD.sizeZ()<<endl;
    QuaternionsOLD.Allocate(RotationsOLD.sizeX(), RotationsOLD.sizeY(),
                            RotationsOLD.sizeZ(), RotationsOLD.Bcells());
//    cout<<"QuaternionsOld allocated yo: Nx = "<<QuaternionsOLD.sizeX()<<"; Ny = "<<QuaternionsOLD.sizeY()<<"; Nz = "<<QuaternionsOLD.sizeZ()<<endl;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, RotationsOLD, RotationsOLD.Bcells(),)
    {
        Quaternion Q;
        Q.set(RotationsOLD(i,j,k));
        QuaternionsOLD(i,j,k) = Q;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    QuaternionsOLD.Remesh(Nx, Ny, Nz);
    RotationsOLD.Reallocate(Nx, Ny, Nz);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, RotationsOLD, RotationsOLD.Bcells(),)
    {
        RotationsOLD(i,j,k) = QuaternionsOLD(i,j,k).RotationMatrix;
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    SetBoundaryConditions(BC);

    Info::WriteStandard(thisclassname, "Remeshed");
}

void LargeDeformations::Solve(vector<OPObject*> Objects, bool screenWrite)
/*
 * Large deformation solver!
 * Calculates the EigenStrain and AppliedStrain increments, and divide them if
 * necessary into small portions, writes the eigenstrain increments into the
 * EP.EffectiveEigenStrain and the AppliedStrain increment into the
 * EP.EffectiveAppliedStrain and calls the elastic solver. The ES.Solve()
 * calculates the corresponding stress and strain increment. From those
 * the velocities are calculated. Using these velocities all the system 
 * variablesis are advected by the deformation inrement. As the second 
 * deformation step, the applied strain increment is considered by the 
 * remeshing and the resulting velocity contribution is added to the velocity 
 * field. Finally the velocity gradient is calculated and the stress is 
 * accumulated considering the local rotations. The whole cycle is repeated 
 * until the total EigenStrain and AppliedStrain increments are applied.
 */
{
    PhaseField* Phase           = (static_cast<PhaseField*>(findOPObject(Objects, "PhaseField", thisclassname, "Solve", true, Verbose)));
    BoundaryConditions* BC      = (static_cast<BoundaryConditions*>(findOPObject(Objects, "BoundaryConditions", thisclassname, "Solve", true, Verbose)));
    ElasticProperties* EP       = (static_cast<ElasticProperties*>(findOPObject(Objects, "ElasticProperties", thisclassname, "Solve", true, Verbose)));
    Velocities* Vel             = (static_cast<Velocities*>(findOPObject(Objects, "Velocities", thisclassname, "Solve", true, Verbose)));
    Orientations* OR            = (static_cast<Orientations*>(findOPObject(Objects, "Orientations", thisclassname, "Solve", true, Verbose)));
    ElasticitySolver* ES        = (static_cast<ElasticitySolver*>(findOPObject(Objects, "ElasticitySolverSpectral", thisclassname, "Solve", true, Verbose)));
    Composition* Cx             = (static_cast<Composition*>(findOPObject(Objects, "Composition", thisclassname, "Solve", false, Verbose)));
    Temperature* Tx             = (static_cast<Temperature*>(findOPObject(Objects, "Temperature", thisclassname, "Solve", false, Verbose)));
    EquilibriumPartitionDiffusion* DF
                                = (static_cast<EquilibriumPartitionDiffusion*>(findOPObject(Objects, "EquilibriumPartitionDiffusion", thisclassname, "Solve", false, Verbose)));
    DrivingForce* dG            = (static_cast<DrivingForce*>(findOPObject(Objects, "DrivingForce", thisclassname, "Solve", false, Verbose)));
//    InterfaceMobility& Mu       = *(static_cast<InterfaceMobility*>(findOPObject(Objects, "InterfaceMobility", thisclassname, "Solve", true, Verbose)));
    InterfaceEnergy* Sigma      = (static_cast<InterfaceEnergy*>(findOPObject(Objects, "InterfaceEnergy", thisclassname, "Solve", true, Verbose)));
    DoubleObstacle* DO            = (static_cast<DoubleObstacle*>(findOPObject(Objects, "DoubleObstacle", thisclassname, "Solve", true, Verbose)));

    DamageModel* DM             = (static_cast<DamageModel*>(findOPObject(Objects, "DamageModel", thisclassname, "Solve", false, Verbose)));
    DamageProperties* DP        = (static_cast<DamageProperties*>(findOPObject(Objects, "DamageProperties", thisclassname, "Solve", false, Verbose)));
    PlasticProperties* PFP    = (static_cast<PlasticProperties*>(findOPObject(Objects, "PlasticProperties", thisclassname, "Solve", false, Verbose)));
    PlasticFlow* PF             = (static_cast<PlasticFlow*>(findOPObject(Objects, "PlasticFlow", thisclassname, "Solve", false, Verbose)));

	int currentIteration = 0;
	int totalSteps = minIterations;
	int newIterations = 1;
	bool force = false;
	if (totalSteps == maxIterations)
	force = true;
	
	/*EP->SetEffectiveEigenStrains(*Phase);
	LargeDeformationsMethods::SaveEffectiveEigenStrains(*Phase, *EP, *OR, *BC,
		EffectiveEigenStrainsOLD);
	EP->AppliedStrainOLD = EP->AppliedStrain;
	EP->AppliedStressOLD = EP->AppliedStress;

	if(DP != nullptr and DM != nullptr)
	{
		EP->SetEffectiveElasticConstants(*Phase);
		EP->ApplyRotation(*OR);
		EP->ApplyDamage(*DP);
	}
	else
	{
		EP->SetEffectiveElasticConstants(*Phase);
		EP->ApplyRotation(*OR);
	}
	LargeDeformationsMethods::SaveEffectiveElasticConstants(*EP, StiffnessOLD);*/

	if (PFP != nullptr) SavePlasticStrain(*PFP, PlasticStrainDamage);

	while (currentIteration < totalSteps)
	{
		if (screenWrite)
			cout << "calculating step " << currentIteration+1 << " of " << totalSteps << " LD step(s)." << endl;

		newIterations = calculateLDStep(Objects, totalSteps, force, screenWrite);
		if (newIterations > 1 and !force)
		{
			if (screenWrite)
			cout << "Step failed, newIt " << newIterations << endl;
			//totalSteps++;// *= 2;
			totalSteps *= 2;
			currentIteration = 0;
			if (totalSteps >= maxIterations)
			{
				force = true;
				totalSteps = maxIterations;
			}
		}
		else
		{
			currentIteration++;
			double LDdt = dt / double(totalSteps);

		//****************************advection:********************************
			bool skip = false;
			if (!skip)
			{
				CalculateAverageVelocity(*EP, *Vel, LDdt);
				Vel->PrescribePhaseVelocities(*Phase);
				Vel->SetBoundaryConditions(VelocityBC);

				double maxVel = Vel->getMaxVelocity();
				int advSubsteps = 1;
				double ADdt = LDdt;
				if (maxVel * LDdt > EP->dx /** 0.5*/)
				{
					advSubsteps = ceil(maxVel*LDdt / (EP->dx * 0.2/** 0.5*/));
					ADdt = LDdt / double(advSubsteps);
					std::string message;
					message = "Advection velocity is too high.\n" +
						Info::GetStandard("Required substeps: ", std::to_string(advSubsteps))
						+ "\n" + Info::GetStandard("Advection subtimestep: ",
							std::to_string(ADdt));
					Info::WriteWarning(message, thisclassname, "Solve()");
				}

				for (int subStep = 0; subStep < advSubsteps; subStep++)
				{
					EP->AdvectStressCompressible(*Vel, *BC, ADdt);
					OR->Advect(*Vel, *BC, ADdt);
					if (Cx != nullptr)
						Cx->AdvectTotal(*Phase, *Vel, *BC, ADdt);
					Phase->Advect(*Vel, *BC, ADdt);
					if (DP != nullptr)
						DP->Advect(*Vel, *BC, ADdt);
					if (PFP != nullptr)
						PFP->Advect(*Vel, *BC, ADdt);
					if (PF != nullptr)
					{
						PlasticFlowCP* tempPF = dynamic_cast<PlasticFlowCP*>(PF);
						if (tempPF != nullptr)
						{
							PlasticFlowMethods::UpdateCRSS(*Phase, *tempPF, *BC);
							tempPF->Advect(*Phase, *Vel, *BC, ADdt);
						}
					}

					Vel->Advect(VelocityBC, ADdt);
					Vel->PrescribePhaseVelocities(*Phase);
					Vel->SetBoundaryConditions(VelocityBC);

					LargeDeformationsMethods::AccumulateStressesRotation(*EP, *OR,
						RotationsOLD, ADdt);
				}
				//**************************end advection*******************************

				//****************************remeshing:********************************
				EP->SetBoundaryConditions(*BC);

				EP->StrainToRemesh += EP->AverageStrain;
				int dN[3];
				dN[0] = EP->StrainToRemesh[0] * Phase->Nx;
				dN[1] = EP->StrainToRemesh[1] * Phase->Ny;
				dN[2] = EP->StrainToRemesh[2] * Phase->Nz;
				int newNx = Phase->Nx + dN[0];
				int newNy = Phase->Ny + dN[1];
				int newNz = Phase->Nz + dN[2];
				if (dN[0] != 0 or dN[1] != 0 or dN[2] != 0)
				{
					stringstream message;
					message << "remeshing..." << endl;
					message << "dN:\t" << dN[0] << "\t" << dN[1] << "\t" << dN[2] << endl;
					Info::WriteSimple(message.str());

					Phase->Remesh(newNx, newNy, newNz, *BC);
					Vel->Remesh(newNx, newNy, newNz, *BC);
					Remesh(newNx, newNy, newNz, *BC);
					EP->Remesh(newNx, newNy, newNz, *BC);
					ES->ReInitialize(*EP);
					OR->Remesh(newNx, newNy, newNz, *BC);
					//RotationsOLD.Reallocate(newNx, newNy, newNz);
					//OR->WriteVTK(0);

					if (Cx != nullptr) Cx->Remesh(newNx, newNy, newNz, *BC);
					if (Tx != nullptr) Tx->Remesh(newNx, newNy, newNz, *BC);
					if (dG != nullptr) dG->Remesh(newNx, newNy, newNz);
					if (DF != nullptr) DF->Remesh(newNx, newNy, newNz);
					if (Sigma != nullptr) Sigma->ReInitialize(*Phase);
					if (DM != nullptr) DM->Reinitialize(newNx, newNy, newNz);
					if (DP != nullptr) DP->Remesh(newNx, newNy, newNz, *BC);
					if (PFP != nullptr) PFP->Remesh(newNx, newNy, newNz, *BC);
					if (PF != nullptr) PF->Remesh(newNx, newNy, newNz, *BC);

					/*vStrain deltaStrain;
					deltaStrain.set_to_zero();
					deltaStrain[0] = dN[0] * 1. / (EP->Nx - dN[0]);
					deltaStrain[1] = dN[1] * 1. / (EP->Ny - dN[1]);
					deltaStrain[2] = dN[2] * 1. / (EP->Nz - dN[2]);
					EP->StrainToRemesh[0] -= deltaStrain[0];
					EP->StrainToRemesh[1] -= deltaStrain[1];
					EP->StrainToRemesh[2] -= deltaStrain[2];*/
				}
				if (IntfStabilisationNeeded and !counter % 10)
				{
					DO->FixSpreading(*Phase, *BC, 10e-6);
				}
				if (PF != nullptr)
				{
					PlasticFlowCP* tempPF = dynamic_cast<PlasticFlowCP*>(PF);
					if (tempPF != nullptr)
					{
						PlasticFlowMethods::UpdateCRSS(*Phase, *tempPF, *BC);
					}
				}

				Vel->PrescribePhaseVelocities(*Phase);
				LargeDeformationsMethods::SetHomogeneousVelocityBC(*Vel);
				Vel->SetBoundaryConditions(VelocityBC);

				LargeDeformationsMethods::CalculateVelocityGradient(*EP, *Vel);
				counter++;
			}
		}
	}
	EP->Energy();
	LargeDeformationsMethods::SaveEffectiveEigenStrains(*Phase, *EP, *OR, *BC,
		EffectiveEigenStrainsOLD);
	EP->AppliedStrainOLD = EP->AppliedStrain;
	EP->AppliedStressOLD = EP->AppliedStress;

	if(DP != nullptr and DM != nullptr)
	{
		EP->SetEffectiveElasticConstants(*Phase);
		EP->ApplyRotation(*OR);
		EP->ApplyDamage(*DP);
	}
	else
	{
		EP->SetEffectiveElasticConstants(*Phase);
		EP->ApplyRotation(*OR);
	}
	LargeDeformationsMethods::SaveEffectiveElasticConstants(*EP, StiffnessOLD);
}

int LargeDeformations::calculateLDStep(vector<OPObject*> Objects, int nSubsteps, bool force, bool screenWrite)
{
	PhaseField* Phase = (static_cast<PhaseField*>(findOPObject(Objects, "PhaseField", thisclassname, "Solve", true, Verbose)));
	BoundaryConditions* BC = (static_cast<BoundaryConditions*>(findOPObject(Objects, "BoundaryConditions", thisclassname, "Solve", true, Verbose)));
	ElasticProperties* EP = (static_cast<ElasticProperties*>(findOPObject(Objects, "ElasticProperties", thisclassname, "Solve", true, Verbose)));
	//Velocities* Vel = (static_cast<Velocities*>(findOPObject(Objects, "Velocities", thisclassname, "Solve", true, Verbose)));
	Orientations* OR = (static_cast<Orientations*>(findOPObject(Objects, "Orientations", thisclassname, "Solve", true, Verbose)));
	ElasticitySolver* ES = (static_cast<ElasticitySolver*>(findOPObject(Objects, "ElasticitySolverSpectral", thisclassname, "Solve", true, Verbose)));
	//Composition* Cx = (static_cast<Composition*>(findOPObject(Objects, "Composition", thisclassname, "Solve", false, Verbose)));
	//Temperature* Tx = (static_cast<Temperature*>(findOPObject(Objects, "Temperature", thisclassname, "Solve", false, Verbose)));
	//EquilibriumPartitionDiffusion* DF
		//= (static_cast<EquilibriumPartitionDiffusion*>(findOPObject(Objects, "EquilibriumPartitionDiffusion", thisclassname, "Solve", false, Verbose)));
	//DrivingForce* dG = (static_cast<DrivingForce*>(findOPObject(Objects, "DrivingForce", thisclassname, "Solve", false, Verbose)));
	//    InterfaceMobility& Mu       = *(static_cast<InterfaceMobility*>(findOPObject(Objects, "InterfaceMobility", thisclassname, "Solve", true, Verbose)));
	//InterfaceEnergy* Sigma = (static_cast<InterfaceEnergy*>(findOPObject(Objects, "InterfaceEnergy", thisclassname, "Solve", true, Verbose)));
	//DoubleObstacle* DO = (static_cast<DoubleObstacle*>(findOPObject(Objects, "DoubleObstacle", thisclassname, "Solve", true, Verbose)));

	DamageModel* DM = (static_cast<DamageModel*>(findOPObject(Objects, "DamageModel", thisclassname, "Solve", false, Verbose)));
	DamageProperties* DP = (static_cast<DamageProperties*>(findOPObject(Objects, "DamageProperties", thisclassname, "Solve", false, Verbose)));
	PlasticProperties* PFP = (static_cast<PlasticProperties*>(findOPObject(Objects, "PlasticProperties", thisclassname, "Solve", false, Verbose)));
	PlasticFlow* PF = (static_cast<PlasticFlow*>(findOPObject(Objects, "PlasticFlow", thisclassname, "Solve", false, Verbose)));

	EP->SetBoundaryConditions(*BC);
	SetBoundaryConditions(*BC);

	double AppliedStrainNorm = LargeDeformationsMethods::GetAppliedStrainNorm(
		*EP, MaxAllowedStrainIncrement,
		screenWrite)
		* nSubsteps;

	double LDdt = dt / double(nSubsteps);

	if (DP != nullptr and DM != nullptr and PFP != nullptr)
	{
		LargeDeformationsMethods::SetPlasticStrainEquivalent(DP->peeq,
			PlasticStrainDamage);

		//            LargeDeformationsMethods::SetPlasticStrainEquivalentIntf(Phase, DP.peeq,
		//                                       CP.AccumulatedEffectiveCreepStrains, 49);
		LargeDeformationsMethods::SaveDamageN(*DP, DamageOLD, DamageNOLD);
		DM->Solve(*EP, *DP, *Phase, *BC);
		OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, DP->EffectiveDamage, 0, )
		{
			double temp = DamageOLD(i, j, k);
			DP->EffectiveDamage(i, j, k) = temp + (DP->EffectiveDamage(i, j, k) - temp)/double(nSubsteps);
			DP->SetBoundaryConditions(*BC);
		}
		OMP_PARALLEL_STORAGE_LOOP_END
		//            DM->DamageStiffness(*EP, *DP);
		EP->SetEffectiveElasticConstants(*Phase);
		EP->ApplyRotation(*OR);
		EP->ApplyDamage(*DP);
	}
	else
	{
		EP->SetEffectiveElasticConstants(*Phase);
		EP->ApplyRotation(*OR);
	}
	LargeDeformationsMethods::SaveRotations(*OR, *BC, RotationsOLD);
	if (PFP != nullptr and PF != nullptr)
	{
		SavePlasticStrain(*PFP, PlasticStrainOLD);
		PlasticFlowCP* tempPF = dynamic_cast<PlasticFlowCP*>(PF);
		if (tempPF != nullptr)
		{
			LargeDeformationsMethods::SaveCRSS(*tempPF, CRSSOLD);
			tempPF->dt = LDdt;
		}
		PF->dt = LDdt;
		PF->Solve(Objects, screenWrite);
		LargeDeformationsMethods::SubtractPlasticStrainOLD(*PFP,
			PlasticStrainOLD);
	}

	double StrainNorm = 0;
	if (!force)
	{
		double EigenStrainNorm = LargeDeformationsMethods::GetEigenstrainNorm(
			Objects, EffectiveEigenStrainsOLD, StiffnessOLD,
			dt, nSubsteps, MaxAllowedStrainIncrement,
			Verbose, screenWrite);

		StrainNorm = min(AppliedStrainNorm, EigenStrainNorm);
	}
	int newIterations = (int)ceil(1. / (StrainNorm));

	if (newIterations == 1 or force)
	{
		LargeDeformationsMethods::SetEigenstrainsIncrement(Objects,
			EffectiveEigenStrainsOLD, StiffnessOLD,
			LDdt, nSubsteps,
			Verbose);

		if (PFP != nullptr and PF != nullptr)
		{
			LargeDeformationsMethods::ADDPlasticStrainOLD(*PFP,
				PlasticStrainOLD);
		}
		LargeDeformationsMethods::SetAppliedStrainIncrement(*EP, nSubsteps);
		LargeDeformationsMethods::SetAppliedStressIncrement(*EP, nSubsteps);
		EP->SetBoundaryConditions(*BC);

		//        int ESiterations =
		ES->Solve(*EP, *OR, *BC, StrainAccuracy, MaxPressure, MaxSolverIterations,
			LDdt);

		return newIterations;
	}
	else
	{
		if (DP != nullptr and DM != nullptr)
		LargeDeformationsMethods::RestoreDamageN(*DP, DamageOLD, DamageNOLD);
		if (PFP != nullptr and PF != nullptr)
		{
			PlasticFlowCP* tempPF = dynamic_cast<PlasticFlowCP*>(PF);
			if (tempPF != nullptr)
			{
				tempPF->dt = LDdt;
				LargeDeformationsMethods::RestoreCRSS(*tempPF, CRSSOLD);
			}
			RestorePlasticStrain(*PFP, PlasticStrainOLD);
		}
		LargeDeformationsMethods::RestoreRotations(*OR, *BC, RotationsOLD);
		return newIterations;
	}
}
//=============================== Private section ==============================

void LargeDeformations::CalculateAverageVelocity(ElasticProperties& EP,
                                                     Velocities& Vel, double dt)
{
    double dtInv = 1.0/dt;
    double dx = EP.dx;
    
    int Nx = EP.Nx;
    int Ny = EP.Ny;
    int Nz = EP.Nz;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EP.Strains, 0,)
    {    
        dVector3 RadiusVector;
        
        RadiusVector[0] = (i - 0.5*Nx)*dx;
        RadiusVector[1] = (j - 0.5*Ny)*dx;       
        RadiusVector[2] = (k - 0.5*Nz)*dx;  
             
        Vel.Average(i,j,k) = (EP.StrainIncrements(i,j,k).tensor() - EP.AverageStrain.tensor())*RadiusVector*dtInv;
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void LargeDeformations::SavePlasticStrain(PlasticProperties& PFP,
                                        Storage3D<vStrain, 0>& PlasticStrainOLD)
{
    int Nx = PFP.Nx;
    int Ny = PFP.Ny;
    int Nz = PFP.Nz;

    if(PlasticStrainOLD.IsNotAllocated())
        PlasticStrainOLD.Allocate(Nx, Ny, Nz, 0);

    if(PlasticStrainOLD.sizeX() != PFP.PlasticStrain.sizeX() or
       PlasticStrainOLD.sizeY() != PFP.PlasticStrain.sizeY() or
       PlasticStrainOLD.sizeZ() != PFP.PlasticStrain.sizeZ())
        PlasticStrainOLD.Reallocate(Nx, Ny, Nz);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, PFP.PlasticStrain, 0,)
    {
        PlasticStrainOLD(i,j,k) = PFP.PlasticStrain(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
void LargeDeformations::RestorePlasticStrain(PlasticProperties& PFP,
                                        Storage3D<vStrain, 0>& PlasticStrainOLD)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, PFP.PlasticStrain, 0,)
    {
        PFP.PlasticStrain(i,j,k) = PlasticStrainOLD(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}
} // namespace opensim
