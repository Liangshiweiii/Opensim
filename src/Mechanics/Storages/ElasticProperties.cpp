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
#include "Mechanics/Storages/ElasticProperties.h"
#include "Mechanics/Storages/SymmetryVariants.h"
#include "Mechanics/Storages/DamageProperties.h"
#include "Mechanics/ElasticityModels/ElasticitySteinbach.h"
#include "Mechanics/ElasticityModels/ElasticityKhachaturyan.h"
#include "Info.h"
#include "PhaseField.h"
#include "Compositions.h"
#include "Orientations.h"
#include "DrivingForce.h"
#include "EquilibriumPartitionDiffusion.h"
#include "Settings.h"
#include "BoundaryConditions.h"
#include "Temperature.h"
#include "Tools/UserInterface.h"
#include "Velocities.h"
#include "VTK.h"
#include "InterfaceEnergy.h"

namespace opensim
{
using namespace std;

ElasticProperties::ElasticProperties(Settings& locSettings,
        const std::string InputFileName)
{
    this->Initialize(locSettings);

    if(InputFileName == "default")
    {
        this->ReadInput();
    }
    else
    {
        this->ReadInput(InputFileName);
    }
}

void ElasticProperties::Initialize(Settings& locSettings)
{
    thisclassname = "ElasticProperties";
    //DefaultInputFileName = ProjectInputDir + "ElasticityInput.opi";

    LargeDeformations = false;
    KeepAspectRatio = false;
    
    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;
    dx = locSettings.dx;
	dt = locSettings.dt;

    ElasticityModel = Khachaturyan;

    Nphases = locSettings.Nphases;
    Ncomp   = max(0,locSettings.Ncomp - 1);
    Names.resize(Ncomp);

    AppStrainMask.set_to_zero();
    AppliedStrain.set_to_zero();
    AppliedStrainBool.set_to_zero();
    AppliedStrainOLD.set_to_zero();
    AppliedStress.set_to_zero();
    AppliedStressOLD.set_to_zero();
    AverageCompliences.set_to_zero();
    AverageElasticConstants.set_to_zero();
    AverageStrain.set_to_zero();
    AvgStrainMask.set_to_zero();
    EffectiveAppliedStrain.set_to_zero();
    EffectiveAppliedStress.set_to_zero();
    LoadStressMask.set_to_zero();
    MAXElasticConstants.set_to_zero();
    MINCompliences.set_to_zero();
    RemeshedStrain.set_to_zero();
    StrainToRemesh.set_to_zero();

    Strains.Allocate(Nx, Ny, Nz, 1);
    StrainIncrements.Allocate(Nx, Ny, Nz, 1);
    Stresses.Allocate(Nx, Ny, Nz, 1);

    EffectiveEigenStrains.Allocate(Nx, Ny, Nz, 1);
    EffectiveElasticConstants.Allocate(Nx, Ny, Nz, 1);

    PhaseElasticConstants.Allocate(Nphases);
    PhaseEigenStrains.Allocate(Nphases);
    PhaseCompliences.Allocate(Nphases);
    PhaseKappa.Allocate({Nphases, Ncomp});
    PhaseLambda.Allocate({Nphases, Ncomp});
    Cref.Allocate({Nphases, Ncomp});

    PhaseAlpha.Allocate(Nphases);
    Tref.Allocate(Nphases);

    ElasticConstants.Allocate(Nphases);
    EigenStrains.Allocate(Nphases);
    Compliences.Allocate(Nphases);

    Kappa.Allocate({Nphases, Ncomp});
    Lambda.Allocate({Nphases, Ncomp});
    Alpha.Allocate(Nphases);

    for(int alpha = 0; alpha != Nphases; alpha++)
    {
        PhaseAlpha[alpha].set_to_zero();
        Tref[alpha] = 0.0;

        PhaseEigenStrains[alpha].set_to_zero();
        PhaseElasticConstants[alpha].set_to_zero();
        PhaseCompliences[alpha].set_to_zero();

        for(int comp = 0; comp < Ncomp; comp++)
        {
            Cref({alpha, comp}) = 0.0;
            PhaseLambda({alpha, comp}).set_to_zero();
            PhaseKappa({alpha, comp}).set_to_zero();
        }
    }

    MAXElasticConstants.set_to_zero();

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Strains,0,)
    {
        Strains(i,j,k).set_to_zero();
        StrainIncrements(i,j,k).set_to_zero();
        Stresses(i,j,k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Variants.Initialize(locSettings);

    Info::WriteStandard("ElasticProperties", "Initialized");
}

void ElasticProperties::InitializeLD()
{
    if (!Stresses.IsAllocated())
    {
        Info::WriteWarning("ElasticProperties should be initialized before LargeDeformations", thisclassname, "InitializedLD");
        exit(1);
    }

    LargeDeformations = true;
    StressIncrements.Allocate(Nx, Ny, Nz, 1);
    VelocityGradient.Allocate(Nx, Ny, Nz, 1);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Strains, 0, )
    {
        StressIncrements(i, j, k).set_to_zero();
        VelocityGradient(i, j, k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Info::WriteStandard("ElasticProperties", "Initialized for Large Deformations");
}

void ElasticProperties::ReadInput(string InputFileName)
{
    Info::WriteLine();
    Info::WriteLineInsert("Elastic properties");
    Info::WriteStandard("Source", InputFileName);

    for(int alpha = 0; alpha != Nphases; alpha++)
    {
        PhaseElasticConstants[alpha].set_to_zero();
    }

    fstream inp(InputFileName.c_str(), ios::in);

    if (!inp)
    {
        Info::WriteExit("File " + InputFileName + " could not be opened",thisclassname, "ReadInput()");
        exit(1);
    };

    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    string tmp1 = UserInterface::ReadParameterF(inp, moduleLocation, "EModel", false, "KHACHATURYAN");
    std::transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::toupper);
    tmp1.erase (std::remove (tmp1.begin(), tmp1.end(), ' '), tmp1.end());

    if (tmp1 == "KHACHATURYAN")
    {
        //Info::WriteStandard("Elasticity model", "Khachaturyan");
        ElasticityModel = Khachaturyan;
    }
    else if (tmp1 == "STEINBACH")
    {
        //Info::WriteStandard("Elasticity model", "Steinbach");
        ElasticityModel = Steinbach;
    }
    else if (tmp1 == "VOIGT")
    {
        //Info::WriteStandard("Elasticity model", "Voigt");
        ElasticityModel = Voigt;
    }
    else if (tmp1 == "REUSS")
    {
        //Info::WriteStandard("Elasticity model", "Reuss");
        ElasticityModel = Reuss;
    }

    /// X direction
    string tmp2 = UserInterface::ReadParameterF(inp, moduleLocation, "BCX", 0);
    std::transform(tmp2.begin(), tmp2.end(), tmp2.begin(), ::toupper);
    tmp2.erase (std::remove (tmp2.begin(), tmp2.end(), ' '), tmp2.end());
    if (tmp2 == "FREEBOUNDARIES")
    {
        //Info::WriteStandard("Boundary Condition X direction", "Free boundary expansion");
        AvgStrainMask[0] = 1.0;
    }
    else if (tmp2 == "APPLIEDSTRAIN")
    {
        //Info::WriteStandard("Boundary Condition X direction", "Applied strain");
        EffectiveAppliedStrain[0] = UserInterface::ReadParameterD(inp, moduleLocation, "BCValueX");
        AppliedStrain[0] = EffectiveAppliedStrain[0];
        AppStrainMask[0] = 1.0;
        Info::WriteStandard("Strain X direction", to_string(EffectiveAppliedStrain[0]));
    }
	else if (tmp2 == "APPLIEDSTRAINRATE")
	{
		//Info::WriteStandard("Boundary Condition X direction", "Strain rate");
		AppliedStrainRate[0] = UserInterface::ReadParameterD(inp, moduleLocation, "BCValueX");
		AppStrainMask[0] = 1.0;
		Info::WriteStandard("Strain rate X direction", to_string(AppliedStrainRate[0]));
	}
    else if (tmp2 == "APPLIEDSTRESS")
    {
        //Info::WriteStandard("Boundary Condition X direction", "Applied stress");
        LoadStressMask[0] = 1.0;
        EffectiveAppliedStress[0] = UserInterface::ReadParameterD(inp, moduleLocation, "BCValueX");
        AppliedStress [0] = EffectiveAppliedStress[0];
        Info::WriteStandard("Stress X direction", to_string(EffectiveAppliedStress[0]));
    }

    /// Y direction
    tmp2 = UserInterface::ReadParameterF(inp, moduleLocation, "BCY", 0);
    std::transform(tmp2.begin(), tmp2.end(), tmp2.begin(), ::toupper);
    tmp2.erase (std::remove (tmp2.begin(), tmp2.end(), ' '), tmp2.end());
    if (tmp2 == "FREEBOUNDARIES")
    {
        //Info::WriteStandard("Boundary Condition Y direction", "Free boundary expansion");
        AvgStrainMask[1] = 1.0;
    }
    else if (tmp2 == "APPLIEDSTRAIN")
    {
        //Info::WriteStandard("Boundary Condition Y direction", "Applied strain");
        EffectiveAppliedStrain[1] = UserInterface::ReadParameterD(inp, moduleLocation, "BCValueY");
        AppliedStrain[1] = EffectiveAppliedStrain[1];
        AppStrainMask[1] = 1.0;
        Info::WriteStandard("Strain Y direction", to_string(EffectiveAppliedStrain[1]));
    }
	else if (tmp2 == "APPLIEDSTRAINRATE")
	{
		//Info::WriteStandard("Boundary Condition X direction", "Strain rate");
		AppliedStrainRate[1] = UserInterface::ReadParameterD(inp, moduleLocation, "BCValueY");
		AppStrainMask[1] = 1.0;
		Info::WriteStandard("Strain rate Y direction", to_string(AppliedStrainRate[1]));
	}
    else if (tmp2 == "APPLIEDSTRESS")
    {
        //Info::WriteStandard("Boundary Condition Y direction", "Applied stress");
        LoadStressMask[1] = 1.0;
        EffectiveAppliedStress[1] = UserInterface::ReadParameterD(inp, moduleLocation, "BCValueY");
        AppliedStress [1] = EffectiveAppliedStress[1];
        Info::WriteStandard("Stress Y direction", to_string(EffectiveAppliedStress[1]));
    }

    /// Z direction
    tmp2 = UserInterface::ReadParameterF(inp, moduleLocation, "BCZ", 0);
    std::transform(tmp2.begin(), tmp2.end(), tmp2.begin(), ::toupper);
    tmp2.erase (std::remove (tmp2.begin(), tmp2.end(), ' '), tmp2.end());
    if (tmp2 == "FREEBOUNDARIES")
    {
        //Info::WriteStandard("Boundary Condition Z direction", "Free boundary expansion");
        AvgStrainMask[2] = 1.0;
    }
    else if (tmp2 == "APPLIEDSTRAIN")
    {
        //Info::WriteStandard("Boundary Condition Z direction", "Applied strain");
        EffectiveAppliedStrain[2] = UserInterface::ReadParameterD(inp, moduleLocation, "BCValueZ");
        AppliedStrain[2] = EffectiveAppliedStrain[2];
        AppStrainMask[2] = 1.0;
        Info::WriteStandard("Strain Z direction", to_string(EffectiveAppliedStrain[2]));
    }
	else if (tmp2 == "APPLIEDSTRAINRATE")
	{
		//Info::WriteStandard("Boundary Condition X direction", "Strain rate");
		AppliedStrainRate[2] = UserInterface::ReadParameterD(inp, moduleLocation, "BCValueZ");
		AppStrainMask[2] = 1.0;
		Info::WriteStandard("Strain rate Y direction", to_string(AppliedStrainRate[2]));
	}
    else if (tmp2 == "APPLIEDSTRESS")
    {
        //Info::WriteStandard("Boundary Condition Z direction", "Applied stress");
        LoadStressMask[2] = 1.0;
        EffectiveAppliedStress[2] = UserInterface::ReadParameterD(inp, moduleLocation, "BCValueZ");
        AppliedStress [2] = EffectiveAppliedStress[2];
        Info::WriteStandard("Stress Z direction", to_string(EffectiveAppliedStress[2]));
    }
    
    string tmp3 = UserInterface::ReadParameterF(inp, moduleLocation, "Restrict", 0);
    std::transform(tmp3.begin(), tmp3.end(), tmp3.begin(), ::toupper);
    tmp3.erase (std::remove (tmp3.begin(), tmp3.end(), ' '), tmp3.end());
    if (tmp3 == "KEEPASPECTRATIO")
    {
        //Info::WriteStandard("Boundary Condition X direction", "Free boundary expansion");
        KeepAspectRatio = true;
    }
    
    // Reading elastic constants
    for(int pIndex = 0; pIndex < Nphases; pIndex++)
    {
        for(int ii =  1; ii <= 6; ii++)
        for(int jj = ii; jj <= 6; jj++)
        {
            stringstream converter;
            converter << "C" << ii << jj << "_" << pIndex;

            PhaseElasticConstants[pIndex](ii-1,jj-1) = UserInterface::ReadParameterD(inp, moduleLocation, converter.str(),false,0);
            if(ii != jj)
            {
                PhaseElasticConstants[pIndex](jj-1,ii-1) = PhaseElasticConstants[pIndex](ii-1,jj-1);
            }
        }
    }
    // Reading eigenstrains
    for(int pIndex = 0; pIndex < Nphases; pIndex++)
    {
        for(int ii =  1; ii <= 6; ii++)
        {
            stringstream converter;
            converter << "E" << ii << "_" << pIndex;

            PhaseEigenStrains[pIndex][ii-1] = UserInterface::ReadParameterD(inp, moduleLocation, converter.str(),false,0);
        }
    }
    // Reading chemo-mechanical coupling parameters if considered
    ChemoMechanicalCoupling = UserInterface::ReadParameterB(inp, moduleLocation, "ChemoMechCoupling");

    if(ChemoMechanicalCoupling)
    {
        for(int comp = 0; comp < Ncomp; comp++)
        {
            stringstream converter;
            converter << string("Comp_") << comp;
            Names[comp] = UserInterface::ReadParameterF(inp, moduleLocation, converter.str());

            for(int alpha = 0; alpha != Nphases; alpha++)
            {
                stringstream converter;
                converter << Names[comp] << "_" << alpha;
                string counter = converter.str();
                Cref({alpha, comp}) = UserInterface::ReadParameterD(inp, moduleLocation, string("CREF_") + counter);
            }
        }

        for (int pIndex = 0; pIndex < Nphases; pIndex++)
        {
            for (int comp = 0; comp < Ncomp; comp++)
            {
                for (int ii = 1; ii <= 6; ii++)
                for (int jj = ii; jj <= 6; jj++)
                {
                    stringstream converter;
                    converter << "Kappa" << ii << jj << "_" << pIndex << "_" << Names[comp];

                    PhaseKappa({ pIndex, comp })(ii - 1, jj - 1) =
                        UserInterface::ReadParameterD(inp, moduleLocation, converter.str(), false, 0.0);
                    if (ii != jj)
                    {
                        PhaseKappa({ pIndex, comp })(jj - 1, ii - 1) =
                            PhaseKappa({ pIndex, comp })(ii - 1, jj - 1);
                    }
                }
            }
        }

        for (int pIndex = 0; pIndex < Nphases; pIndex++)
        {
            for (int comp = 0; comp < Ncomp; comp++)
            {
                for (int ii = 1; ii <= 6; ii++)
                {
                    stringstream converter;
                    converter << "Lambda" << ii << "_" << pIndex << "_" << Names[comp];

                    PhaseLambda({ pIndex, comp })[ii - 1] =
                        UserInterface::ReadParameterD(inp, moduleLocation, converter.str(), false, 0.0);
                }
            }
        }
    }
    // Reading thermo-mechanical coupling parameters if considered
    ThermoMechanicalCoupling = UserInterface::ReadParameterB(inp, moduleLocation, "ThermoMechCoupling");
    if(ThermoMechanicalCoupling)
    {
        //TODO: implement reading of thermo-mechanical coupling parameters
    }

	NeuberCorrection = UserInterface::ReadParameterB(inp, moduleLocation, "NeuberCorrection", false, "No");
	if (NeuberCorrection)
	{
		//TODO: implement reading of thermo-mechanical coupling parameters
	}
    inp.close();

    stringstream outstream;

    Info::WriteBlankLine();
    for (int pIndex = 0; pIndex < Nphases; pIndex++)
    {
        outstream << "Elastic properties of Phase " << pIndex << ":" << endl;
        // Check if stiffness constants are read.
        if (PhaseElasticConstants[pIndex].norm() == 0)
        {
            std::string message = "Elastic constants for phase " + to_string(pIndex) + " not set.";
            Info::WriteExit(message, thisclassname, "ReadInput()");
            exit(3);
        }
        else
        {
            PhaseCompliences[pIndex] = PhaseElasticConstants[pIndex].inverted();
        }
        outstream << "E_" << pIndex << ": " << PhaseEigenStrains[pIndex].print() << endl;
        outstream << "C_" << pIndex << ": " << endl << PhaseElasticConstants[pIndex].print() << endl;
        for(int comp = 0; comp < Ncomp; comp++)
        {
            if (PhaseKappa({pIndex, comp}).norm() != 0.0)
            {
                outstream << "Kappa_" << pIndex << "_" << Names[comp] << ": " << endl << PhaseKappa({pIndex, comp}).print() << endl;
            }
            if (PhaseLambda({pIndex, comp}).norm() != 0.0)
            {
                outstream << "Lambda_" << pIndex << "_" << Names[comp] << ": " << PhaseLambda({pIndex, comp}).print() << endl;
            }
        }
        if (PhaseAlpha[pIndex].norm() != 0.0)
        {
            outstream << "Alpha_" << pIndex << ": " << PhaseAlpha[pIndex].print() << endl;
            outstream << "Tref: " << Tref[pIndex] << endl;
        }
    }
    Info::WriteSimple(outstream.str());
    Info::WriteLine();
}

void ElasticProperties::Remesh(int newNx, int newNy, int newNz, BoundaryConditions& BC)
{
    ///Changing box size
    double strainInc[3];
    strainInc[0] = double(newNx - Nx)/Nx;
    strainInc[1] = double(newNy - Ny)/Ny;
    strainInc[2] = double(newNz - Nz)/Nz;

    RemeshedStrain[0] += strainInc[0];
    RemeshedStrain[1] += strainInc[1];
    RemeshedStrain[2] += strainInc[2];

    StrainToRemesh[0] -= strainInc[0];
    StrainToRemesh[1] -= strainInc[1];
    StrainToRemesh[2] -= strainInc[2];

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    Stresses.Remesh(Nx, Ny, Nz);
    StressIncrements.Remesh(Nx, Ny, Nz);
    Strains.Remesh(Nx, Ny, Nz);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, Strains, 0,)
    {
        Strains(i,j,k)[0] += strainInc[0];
        Strains(i,j,k)[1] += strainInc[1];
        Strains(i,j,k)[2] += strainInc[2];
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    StrainIncrements.Remesh(Nx, Ny, Nz);
//    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, Strains, 0,)
//    {
//        StrainIncrements(i,j,k)[0] += strainInc[0];
//        StrainIncrements(i,j,k)[1] += strainInc[1];
//        StrainIncrements(i,j,k)[2] += strainInc[2];
//    }
//    OMP_PARALLEL_STORAGE_LOOP_END
//    AverageStrain[0] += strainInc[0];
//    AverageStrain[1] += strainInc[1];
//    AverageStrain[2] += strainInc[2];

    EffectiveEigenStrains.Remesh(Nx, Ny, Nz);
    EffectiveElasticConstants.Remesh(Nx, Ny, Nz);

    if(VelocityGradient.IsAllocated()) VelocityGradient.Reallocate(Nx, Ny, Nz);

    SetBoundaryConditions(BC);

    //WriteRemeshingData(tStep, sim_time);
    Info::WriteStandard(thisclassname, "Remeshed");
}

void ElasticProperties::SetGrainsProperties(PhaseField& Phase)
{
    int size = Phase.FieldsStatistics.size();
    if(EigenStrains.Size() != size)
    {
        EigenStrains.Reallocate(size);
        ElasticConstants.Reallocate(size);
        Compliences.Reallocate(size);

        Lambda.Reallocate({size, Ncomp});
        Kappa.Reallocate({size, Ncomp});
        Alpha.Reallocate(size);
    }

    for(int alpha = 0; alpha != size; alpha++)
    if(Phase.FieldsStatistics[alpha].Exist)
    {
        int pIndex = Phase.FieldsStatistics[alpha].Phase;
        int vIndex = Phase.FieldsStatistics[alpha].Variant;

        EigenStrains[alpha]     = PhaseEigenStrains[pIndex];
        ElasticConstants[alpha] = PhaseElasticConstants[pIndex];
        Alpha[alpha] = PhaseAlpha[pIndex];

        if(Variants.set)
        {
            EigenStrains[alpha].rotate(Variants(pIndex, vIndex));
            ElasticConstants[alpha].rotate(Variants(pIndex, vIndex));
            Alpha[alpha].rotate(Variants(pIndex, vIndex));
        }

        EigenStrains[alpha].rotate(Phase.FieldsStatistics[alpha].Orientation.RotationMatrix);
        ElasticConstants[alpha].rotate(Phase.FieldsStatistics[alpha].Orientation.RotationMatrix);
        Alpha[alpha].rotate(Phase.FieldsStatistics[alpha].Orientation.RotationMatrix);

        Compliences[alpha] = ElasticConstants[alpha].inverted();

        for(int comp = 0; comp  < Ncomp; comp++)
        {
            Lambda({alpha, comp}) = PhaseLambda({pIndex, comp});
            Kappa({alpha, comp})  = PhaseKappa({pIndex, comp});

            if(Variants.set)
            {
                Lambda({alpha, comp}).rotate(Variants(pIndex, vIndex));
                Kappa({alpha, comp}).rotate(Variants(pIndex, vIndex));
            }

            Lambda({alpha, comp}).rotate(Phase.FieldsStatistics[alpha].Orientation.RotationMatrix);
            Kappa({alpha, comp}).rotate(Phase.FieldsStatistics[alpha].Orientation.RotationMatrix);
        }
    }

    for(int alpha = 0; alpha != size; alpha++)
    if(Phase.FieldsStatistics[alpha].Exist)
    for(int n = 0; n < 6; n++)
    for(int m = 0; m < 6; m++)
    {
        MAXElasticConstants(n,m) = std::max(ElasticConstants[alpha](n,m), MAXElasticConstants(n,m));
        MINCompliences(n,m) = std::min(Compliences[alpha](n,m), MINCompliences(n,m));
    }
}

void ElasticProperties::SetBoundaryConditions(const BoundaryConditions& BC)
{
    // Only periodic BC are correct. For non periodic boundary conditions gradients should be treated differently.

    BC.SetX(Strains);
    BC.SetY(Strains);
    BC.SetZ(Strains);

    BC.SetX(Stresses);
    BC.SetY(Stresses);
    BC.SetZ(Stresses);

    BC.SetX(StrainIncrements);
    BC.SetY(StrainIncrements);
    BC.SetZ(StrainIncrements);

    BC.SetX(EffectiveEigenStrains);
    BC.SetY(EffectiveEigenStrains);
    BC.SetZ(EffectiveEigenStrains);

    BC.SetX(EffectiveElasticConstants);
    BC.SetY(EffectiveElasticConstants);
    BC.SetZ(EffectiveElasticConstants);

    if(LargeDeformations)
    {
        BC.SetX(StressIncrements);
        BC.SetY(StressIncrements);
        BC.SetZ(StressIncrements);

        BC.SetX(VelocityGradient);
        BC.SetY(VelocityGradient);
        BC.SetZ(VelocityGradient);
    }
}

ElasticProperties& ElasticProperties::operator= (const ElasticProperties& rhs)
{
    // protect against invalid self-assignment and copy of unitialized object
    if (this != &rhs and rhs.thisclassname == "ElasticProperties")
    {
        thisclassname = rhs.thisclassname;
        //DefaultInputFileName = rhs.DefaultInputFileName;

        LargeDeformations = rhs.LargeDeformations;

        Nx = rhs.Nx;
        Ny = rhs.Ny;
        Nz = rhs.Nz;
        dx = rhs.dx;
		dt = rhs.dt;

        Nphases = rhs.Nphases;
        Ncomp   = rhs.Ncomp;
        Names = rhs.Names;
        ElasticityModel = rhs.ElasticityModel;

        if (Strains.IsNotAllocated())
        {
            Strains.Allocate(Nx, Ny, Nz, rhs.Strains.Bcells());
            StrainIncrements.Allocate(Nx, Ny, Nz, rhs.StrainIncrements.Bcells());
            Stresses.Allocate(Nx, Ny, Nz, rhs.Stresses.Bcells());
            EffectiveEigenStrains.Allocate(Nx, Ny, Nz, rhs.EffectiveEigenStrains.Bcells());
            EffectiveElasticConstants.Allocate(Nx, Ny, Nz, rhs.EffectiveElasticConstants.Bcells());

            if (LargeDeformations)
            {
                VelocityGradient.Allocate(Nx, Ny, Nz, rhs.VelocityGradient.Bcells());
                StressIncrements.Allocate(Nx, Ny, Nz, rhs.StressIncrements.Bcells());
            }
        }
        else if (not Strains.IsSize(rhs.Nx, rhs.Ny, rhs.Nz))
        {
            Strains.Reallocate(Nx, Ny, Nz);
            StrainIncrements.Reallocate(Nx, Ny, Nz);
            Stresses.Reallocate(Nx, Ny, Nz);
            EffectiveEigenStrains.Reallocate(Nx, Ny, Nz);
            EffectiveElasticConstants.Reallocate(Nx, Ny, Nz);

            if (LargeDeformations)
            {
                StressIncrements.Reallocate(Nx, Ny, Nz);
                VelocityGradient.Reallocate(Nx, Ny, Nz);
            }
        }

        if (PhaseElasticConstants.IsNotAllocated())
        {
            PhaseElasticConstants.Allocate(Nphases);
            PhaseEigenStrains.Allocate(Nphases);
            PhaseCompliences.Allocate(Nphases);
            PhaseKappa.Allocate({Nphases, Ncomp});
            PhaseLambda.Allocate({Nphases, Ncomp});
            PhaseAlpha.Allocate(Nphases);
            Tref.Allocate(Nphases);
            Cref.Allocate({Nphases, Ncomp});
        }
        else if (PhaseElasticConstants.Size() != rhs.Nphases)
        {
            PhaseElasticConstants.Reallocate(Nphases);
            PhaseEigenStrains.Reallocate(Nphases);
            PhaseCompliences.Reallocate(Nphases);
            PhaseKappa.Reallocate({Nphases, Ncomp});
            PhaseLambda.Reallocate({Nphases, Ncomp});
            PhaseAlpha.Reallocate(Nphases);
            Tref.Reallocate(Nphases);
            Cref.Reallocate({Nphases, Ncomp});
        }
        for(int alpha = 0; alpha != Nphases; alpha++)
        {
            PhaseAlpha[alpha] = rhs.PhaseAlpha[alpha];
            Tref[alpha] = rhs.Tref[alpha];

            PhaseEigenStrains[alpha] = rhs.PhaseEigenStrains[alpha];
            PhaseElasticConstants[alpha] = rhs.PhaseElasticConstants[alpha];
            PhaseCompliences[alpha] = rhs.PhaseCompliences[alpha];

            for(int comp = 0; comp < Ncomp; comp++)
            {
                PhaseLambda({alpha, comp}) = rhs.PhaseLambda({alpha, comp});
                PhaseKappa({alpha, comp}) = rhs.PhaseKappa({alpha, comp});
                Cref({alpha, comp}) = rhs.Cref({alpha, comp});
            }
        }
        int n = 0;
        if (rhs.ElasticConstants.Size() > rhs.Nphases) // if SetGrainsProperties has been called
        {
            n = rhs.ElasticConstants.Size();
        }
        else
        {
            n = rhs.Nphases;
        }
        if (ElasticConstants.IsNotAllocated())
        {
            ElasticConstants.Allocate(n);
            EigenStrains.Allocate(n);
            Compliences.Allocate(n);
            Kappa.Allocate({n, Ncomp});
            Lambda.Allocate({n, Ncomp});
            Alpha.Allocate(n);
        }
        else if (ElasticConstants.Size() != n)
        {
            EigenStrains.Reallocate(n);
            ElasticConstants.Reallocate(n);
            Compliences.Reallocate(n);
            Lambda.Reallocate({n, Ncomp});
            Kappa.Reallocate({n, Ncomp});
            Alpha.Reallocate(n);
        }
        for(int alpha = 0; alpha != n; alpha++)
        {
            ElasticConstants[alpha] = rhs.ElasticConstants[alpha];
            EigenStrains[alpha] = rhs.EigenStrains[alpha];
            Compliences[alpha] = rhs.Compliences[alpha];
            Alpha[alpha] = rhs.Alpha[alpha];
            for(int comp = 0; comp < Ncomp; comp++)
            {
                Lambda({alpha, comp}) = rhs.Lambda({alpha, comp});
                Kappa({alpha, comp}) = rhs.Kappa({alpha, comp});
            }
        }

        RemeshedStrain = rhs.RemeshedStrain;
        StrainToRemesh = rhs.StrainToRemesh;
        AverageStrain = rhs.AverageStrain;
        AppliedStrainBool = rhs.AppliedStrainBool;
        AppliedStress = rhs.AppliedStress;
        AppliedStrain = rhs.AppliedStrain;
        AppliedStrainOLD = rhs.AppliedStrainOLD;
        EffectiveAppliedStrain = rhs.EffectiveAppliedStrain;
        EffectiveAppliedStress = rhs.EffectiveAppliedStress;
        MAXElasticConstants = rhs.MAXElasticConstants;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Strains,Strains.Bcells(),)
        {
            Strains(i,j,k) = rhs.Strains(i,j,k);
            StrainIncrements(i,j,k) = rhs.StrainIncrements(i,j,k);
            Stresses(i,j,k) = rhs.Stresses(i,j,k);
            EffectiveEigenStrains(i,j,k) = rhs.EffectiveEigenStrains(i,j,k);
            EffectiveElasticConstants(i,j,k) = rhs.EffectiveElasticConstants(i,j,k);
            if (LargeDeformations)
            {
                StressIncrements(i, j, k) = rhs.StressIncrements(i, j, k);
                VelocityGradient(i, j, k) = rhs.VelocityGradient(i, j, k);
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    return *this;
}

void ElasticProperties::CalculateThermalStrains(PhaseField& Phase, Temperature& Tx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        {
            int pIndex = Phase.FieldsStatistics[alpha->index].Phase;
            for(int n = 0; n < 3; n++)
            {
                EffectiveEigenStrains(i,j,k)[n] +=
                        Alpha[alpha->index][n] * (Tx.Tx(i,j,k) - Tref[pIndex])*(alpha->value);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}


void ElasticProperties::SetEffectiveEigenStrains(PhaseField& Phase)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EffectiveEigenStrains, EffectiveEigenStrains.Bcells(),)
    {
        EffectiveEigenStrains(i,j,k).set_to_zero();
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha < Phase.Fields(i,j,k).cend(); ++alpha)
        {
            EffectiveEigenStrains(i,j,k) += EigenStrains[alpha->index]*alpha->value;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::SetEffectiveEigenStrains(PhaseField& Phase, Composition& Cx)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k, EffectiveEigenStrains, EffectiveEigenStrains.Bcells(),)
    {
        EffectiveEigenStrains(i,j,k).set_to_zero();

        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha < Phase.Fields(i,j,k).cend(); ++alpha)
        {
            int index = alpha->index;
            int pIndex = Phase.FieldsStatistics[index].Phase;
            double phiAlpha = alpha->value;

            EffectiveEigenStrains(i,j,k) += EigenStrains[index]*phiAlpha;

            for(int comp = 0; comp < Ncomp; comp++)
            {
                double delta = (Cx.Phase(i,j,k)({pIndex, comp}) -
                                               Cref({pIndex, comp}));

                for(int dir = 0; dir < 6; dir++)
                {
                    EffectiveEigenStrains(i,j,k)[dir] += delta*
                                 phiAlpha*Lambda({index, comp})[dir];
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void ElasticProperties::SetEffectiveElasticConstants(PhaseField& Phase)
{
    switch(ElasticityModel)
    {
        case Khachaturyan:
        {
            ElasticityKhachaturyan::SetEffectiveElasticConstants(Phase, *this);
            break;
        }
        case Steinbach:
        {
            ElasticitySteinbach::SetEffectiveElasticConstants(Phase, *this);
            break;
        }
        default:
        {
            Info::WriteExit("Nonexistent elasticity model selected", thisclassname, "SetEffectiveElasticConstants()");
            exit(1);
        }
    }
}

void ElasticProperties::SetEffectiveElasticConstants(PhaseField& Phase, Composition& Cx)
{
    switch(ElasticityModel)
    {
        case Khachaturyan:
        {
            ElasticityKhachaturyan::SetEffectiveElasticConstants(Phase, *this, Cx);
            break;
        }
        case Steinbach:
        {
            ElasticitySteinbach::SetEffectiveElasticConstants(Phase, *this, Cx);
            break;
        }
        default:
        {
            Info::WriteExit("Nonexistent elasticity model selected", thisclassname, "SetEffectiveElasticConstants()");
            exit(1);
        }
    }
}

void ElasticProperties::CalculateDrivingForce(PhaseField& Phase, DrivingForce& dGab)
{
    switch(ElasticityModel)
    {
        case Khachaturyan:
        {
            ElasticityKhachaturyan::CalculateDrivingForce(Phase, *this, dGab);
            break;
        }
        case Steinbach:
        {
            ElasticitySteinbach::CalculateDrivingForce(Phase, *this, dGab);
            break;
        }
        default:
        {
            Info::WriteExit("Nonexistent elasticity model selected", thisclassname, "CalculateDrivingForce()");
            exit(1);
        }
    }
}

void ElasticProperties::CalculateChemicalPotentialContribution(PhaseField& Phase,
                                              EquilibriumPartitionDiffusion& DF)
{
    switch(ElasticityModel)
    {
        case Khachaturyan:
        {
            ElasticityKhachaturyan::CalculateChemicalPotentialContribution(Phase, *this, DF);
            break;
        }
        case Steinbach:
        {
            ElasticitySteinbach::CalculateChemicalPotentialContribution(Phase, *this, DF);
            break;
        }
        default:
        {
            Info::WriteExit("Nonexistent elasticity model selected", thisclassname, "CalculateChemicalPotentialContribution()");
            exit(1);
        }
    }
}

void ElasticProperties::AdvectEigenStrain(Velocities& Vel, BoundaryConditions& BC, double dt, int scheme)
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
                     ((fabs(Vel.Average(i-1,j,k)[0]) + Vel.Average(i-1,j,k)[0])*Strains(i-1,j,k)[n] +
                      (fabs(Vel.Average(i,j-1,k)[1]) + Vel.Average(i,j-1,k)[1])*Strains(i,j-1,k)[n] +
                      (fabs(Vel.Average(i,j,k-1)[2]) + Vel.Average(i,j,k-1)[2])*Strains(i,j,k-1)[n] +
                      (fabs(Vel.Average(i+1,j,k)[0]) - Vel.Average(i+1,j,k)[0])*Strains(i+1,j,k)[n] +
                      (fabs(Vel.Average(i,j+1,k)[1]) - Vel.Average(i,j+1,k)[1])*Strains(i,j+1,k)[n] +
                      (fabs(Vel.Average(i,j,k+1)[2]) - Vel.Average(i,j,k+1)[2])*Strains(i,j,k+1)[n]) -
                      (fabs(Vel.Average(i,j,k)[0]) +
                       fabs(Vel.Average(i,j,k)[1]) +
                       fabs(Vel.Average(i,j,k)[2])) * Strains(i, j, k)[n]/dx;
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            break;
        }
        case LaxWendroff:
        {
            Info::WriteExit("LaxWendroff advection scheme is not supported for EigenStrains",
                            thisclassname, "AdvectEigenStrain(EP, Vel, BC, dt)");
            exit(13);
            break;
        }
        default:
        {
            Info::WriteExit("Wrong/No advection scheme is given in the input file",
                             thisclassname, "AdvectEigenStrain(EP, Vel, BC, dt)");
            exit(13);
        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,StrainsDot,0,)
    {
        for (int n = 0; n < 6; ++n)
        {
            EffectiveEigenStrains(i, j, k)[n] += StrainsDot(i, j, k)[n]*dt;
            StrainsDot(i, j, k)[n] = 0.0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    SetBoundaryConditions(BC);
}

void ElasticProperties::AdvectStress(Velocities& Vel, BoundaryConditions& BC, double dt, int scheme)
{
    if(StressesDot.IsNotAllocated())
    {
        StressesDot.Allocate(Nx, Ny, Nz, 1);
    }
    if(not StressesDot.IsSize(Nx, Ny, Nz))
    {
        StressesDot.Reallocate(Nx, Ny, Nz);
    }
    SetBoundaryConditions(BC);

    switch(scheme)
    {
        case Upwind:
        {
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,StressesDot,0,)
            for (int n = 0; n < 6; n++)
            {
                StressesDot(i,j,k)[n] =
                    -Vel.Average(i,j,k)[0]/dx*
                    ((Stresses(i+1, j, k)[n] - Stresses(i, j, k)[n])*(Vel.Average(i,j,k)[0]<0.0) +
                    (-Stresses(i-1, j, k)[n] + Stresses(i, j, k)[n])*(Vel.Average(i,j,k)[0]>0.0)) -
                    Vel.Average(i,j,k)[1]/dx*
                    ((Stresses(i, j+1, k)[n] - Stresses(i, j, k)[n])*(Vel.Average(i,j,k)[1]<0.0) +
                    (-Stresses(i, j-1, k)[n] + Stresses(i, j, k)[n])*(Vel.Average(i,j,k)[1]>0.0)) -
                    Vel.Average(i,j,k)[2]/dx*
                    ((Stresses(i, j, k+1)[n] - Stresses(i, j, k)[n])*(Vel.Average(i,j,k)[2]<0.0) +
                    (-Stresses(i, j, k-1)[n] + Stresses(i, j, k)[n])*(Vel.Average(i,j,k)[2]>0.0));
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            break;
        }
        case LaxWendroff:
        {
            Info::WriteExit("LaxWendroff advection scheme is not supported for AccStress",
                            thisclassname, "AdvectAccStress(EP, Vel, BC, dt)");
            exit(13);
            break;
        }
        default:
        {
            Info::WriteExit("Wrong/No advection scheme is given in the input file",
                             thisclassname, "AdvectAccStress(EP, Vel, BC, dt)");
            exit(13);
        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,StressesDot,0,)
    {
        for (int n = 0; n < 6; ++n)
        {
            Stresses(i, j, k)[n] += StressesDot(i, j, k)[n]*dt;
            StressesDot(i, j, k)[n] = 0.0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    SetBoundaryConditions(BC);
}

void ElasticProperties::AdvectStressCompressible(Velocities& Vel, BoundaryConditions& BC, double dt, int scheme)
{
    if(StressesDot.IsNotAllocated())
    {
        StressesDot.Allocate(Nx, Ny, Nz, 1);
    }
    if(not StressesDot.IsSize(Nx, Ny, Nz))
    {
        StressesDot.Reallocate(Nx, Ny, Nz);
    }
    SetBoundaryConditions(BC);

    switch(scheme)
    {
    case Upwind:
    {
        const double dx2Inv = 0.5/dx;
        const double dxInv = 1.0/dx;
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,StressesDot,0,)
        for (int n = 0; n < 6; n++)
        {
            double val   = Stresses(i,j,k)[n];
            double valXm = Stresses(i-1,j,k)[n];
            double valXp = Stresses(i+1,j,k)[n];
            double valYm = Stresses(i,j-1,k)[n];
            double valYp = Stresses(i,j+1,k)[n];
            double valZm = Stresses(i,j,k-1)[n];
            double valZp = Stresses(i,j,k+1)[n];

            double ux = Vel.Average(i,j,k)[0];
            double uy = Vel.Average(i,j,k)[1];
            double uz = Vel.Average(i,j,k)[2];

            double uxm = Vel.Average(i-1,j,k)[0];
            double uxp = Vel.Average(i+1,j,k)[0];

            double uym = Vel.Average(i,j-1,k)[1];
            double uyp = Vel.Average(i,j+1,k)[1];

            double uzm = Vel.Average(i,j,k-1)[2];
            double uzp = Vel.Average(i,j,k+1)[2];

            StressesDot(i,j,k)[n] = dx2Inv*
                    ((fabs(uxm) + uxm)*valXm +
                     (fabs(uym) + uym)*valYm +
                     (fabs(uzm) + uzm)*valZm +
                     (fabs(uxp) - uxp)*valXp +
                     (fabs(uyp) - uyp)*valYp +
                     (fabs(uzp) - uzp)*valZp) -
                     (fabs(ux) + fabs(uy) + fabs(uz)) * val * dxInv;
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        break;
    }
    case LaxWendroff:
    {
        Info::WriteExit("LaxWendroff advection scheme is not supported for AccStress",
                        thisclassname, "AdvectAccStress(EP, Vel, BC, dt)");
        exit(13);
        break;
    }
    default:
    {
        Info::WriteExit("Wrong/No advection scheme is given in the input file",
                         thisclassname, "AdvectAccStress(EP, Vel, BC, dt)");
        exit(13);
    }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,StressesDot,0,)
    {
        for (int n = 0; n < 6; ++n)
        {
            Stresses(i, j, k)[n] += StressesDot(i, j, k)[n]*dt;
            StressesDot(i, j, k)[n] = 0.0;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    SetBoundaryConditions(BC);
}

void ElasticProperties::WriteStrains(int tStep)
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"Strains_", tStep, ".dat");

    fstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        Info::WriteExit("File " + FileName + " could not be created",thisclassname);
        exit(1);
    };

    STORAGE_LOOP_BEGIN(i,j,k,Strains,Strains.Bcells())
    {
        out.write(reinterpret_cast<char*>(&(Strains(i,j,k))), sizeof(vStrain));
    }
    STORAGE_LOOP_END

    out.close();
}

void ElasticProperties::WriteStresses(int tStep)
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"Stresses_", tStep, ".dat");

    fstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        Info::WriteExit("File " + FileName + " could not be created",thisclassname);
        exit(1);
    };

    STORAGE_LOOP_BEGIN(i,j,k,Stresses,Stresses.Bcells())
    {
        out.write(reinterpret_cast<char*>(&(Stresses(i,j,k))), sizeof(vStress));
    }
    STORAGE_LOOP_END
    out.close();
}

void ElasticProperties::ReadStrains(int tStep)
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"Strains_", tStep, ".dat");

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit("File " + FileName + " could not be opened",thisclassname);
        exit(1);
    };

    STORAGE_LOOP_BEGIN(i,j,k,Strains,Strains.Bcells())
    {
        inp.read(reinterpret_cast<char*>(&(Strains(i,j,k))), sizeof(vStrain));
    }
    STORAGE_LOOP_END
    Info::WriteStandard("Strains", "Binary input read successfully");
}

void ElasticProperties::ReadStresses(int tStep)
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"Stresses_", tStep, ".dat");

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit("File " + FileName + " could not be opened",thisclassname);
        exit(1);
    };

    STORAGE_LOOP_BEGIN(i,j,k,Stresses,Stresses.Bcells())
    {
        inp.read(reinterpret_cast<char*>(&(Stresses(i,j,k))), sizeof(vStress));
    }
    STORAGE_LOOP_END
}

void ElasticProperties::WriteStrainsVTK(const int tStep) const
{
    std::stringstream buffer;
    std::vector<int> DataTypes {PDScalars};
    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        WriteStrainsVTKData(buffer);
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "Strains", tStep);
}

void ElasticProperties::WriteStrainsVTKData(std::stringstream& buffer) const
{
    vector<int> compV {0, 1, 2, 3, 4, 5};
    vector<string> compNameV {"1", "2", "3", "4", "5", "6"};

    for (auto it = compV.cbegin(); it != compV.cend(); ++it)
    {
        buffer << "<DataArray type = \"Float64\" Name = \""
               << "E_" << compNameV[*it]
               << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            buffer << Strains(i,j,k)[*it] << endl;
        }
        buffer << "</DataArray>" << endl;
    }
}

void ElasticProperties::WriteEffectiveEigenStrainsTensorVTK(int tStep)
{
    int thisStorageSizeX = Strains.sizeX();
    int thisStorageSizeY = Strains.sizeY();
    int thisStorageSizeZ = Strains.sizeZ();

    stringstream buffer;
    std::vector<int> DataTypes {PDScalars};
    vector<long int> dimV {thisStorageSizeX, thisStorageSizeY, thisStorageSizeZ};
    vector<int> compV {0, 1, 2, 3, 4, 5};
    vector<string> compNameV {"xx", "yy", "zz", "yz", "xz", "xy"};

    VTK::WriteHeader(buffer, thisStorageSizeX, thisStorageSizeY, thisStorageSizeZ);
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        for (auto it = compV.cbegin(); it != compV.cend(); ++it)
        {
            if(*it < 3)
            {
                buffer << "<DataArray type = \"Float64\" Name = \""
                       << "EpsilonEigen_" << compNameV[*it]
                       << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
                for(int k = 0; k < thisStorageSizeZ; ++k)
                for(int j = 0; j < thisStorageSizeY; ++j)
                for(int i = 0; i < thisStorageSizeX; ++i)
                {
                    buffer << EffectiveEigenStrains(i,j,k)[*it] << endl;
                }
                buffer << "</DataArray>" << endl;
            }
            else
            {
                buffer << "<DataArray type = \"Float64\" Name = \""
                       << "EpsilonEigen_" << compNameV[*it]
                       << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
                for(int k = 0; k < thisStorageSizeZ; ++k)
                for(int j = 0; j < thisStorageSizeY; ++j)
                for(int i = 0; i < thisStorageSizeX; ++i)
                {
                    // Remove factor 2 of Voigt notation
                    buffer << 0.5*EffectiveEigenStrains(i,j,k)[*it] << endl;
                }
                buffer << "</DataArray>" << endl;
            }
        }
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, thisStorageSizeX, thisStorageSizeY, thisStorageSizeZ);
    VTK::WriteToFile(buffer, "EffEigenStrains", tStep);
}

void ElasticProperties::WriteEffectiveEigenStrainsVoigtVTK(int tStep)
{
    stringstream buffer;
    std::vector<int> DataTypes {PDScalars};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    for(int m = 0; m < 6; m++)
    {
        buffer << "<DataArray type = \"Float64\" Name = \""
               << "E_" << m + 1
               << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < EffectiveEigenStrains.sizeZ(); ++k)
        for(int j = 0; j < EffectiveEigenStrains.sizeY(); ++j)
        for(int i = 0; i < EffectiveEigenStrains.sizeX(); ++i)
        {
            buffer << EffectiveEigenStrains(i,j,k)[m] << endl;
        }
        buffer << "</DataArray>" << endl;
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "EffectiveEigenStrains", tStep);
}

void ElasticProperties::WriteEffectiveElasticConstantsVTK(int tStep)
{
    stringstream buffer;
    std::vector<int> DataTypes {PDScalars};
    vector<long int> dimV {EffectiveEigenStrains.sizeX(),
                  EffectiveEigenStrains.sizeY(), EffectiveEigenStrains.sizeZ()};

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
            buffer << EffectiveElasticConstants(i,j,k)(m,n) << endl;
        }
        buffer << "</DataArray>" << endl;
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "EffectiveElasticConstants", tStep);
}

void ElasticProperties::WriteStressesVTK(const int tStep) const
{
    std::stringstream buffer;
    std::vector<int> DataTypes {PDScalars};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        WriteStressesVTKData(buffer);
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "Stresses", tStep);
}

void ElasticProperties::WriteStressesVTKData(std::stringstream& buffer) const
{
    VTK::WriteSymmetricTensor(buffer, Stresses, "Sigma");

    buffer << "<DataArray type = \"Float64\" Name = \"" << "Pressure" <<
                "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for(int k = 0; k < Stresses.sizeZ(); ++k)
    for(int j = 0; j < Stresses.sizeY(); ++j)
    for(int i = 0; i < Stresses.sizeX(); ++i)
    {
        buffer << Stresses(i,j,k).Pressure() << endl;
    }
    buffer << "</DataArray>" << endl;

    buffer << "<DataArray type = \"Float64\" Name = \"" << "vMises" <<
                "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for(int k = 0; k < Stresses.sizeZ(); ++k)
    for(int j = 0; j < Stresses.sizeY(); ++j)
    for(int i = 0; i < Stresses.sizeX(); ++i)
    {
        buffer << Stresses(i,j,k).Mises() << endl;
    }
    buffer << "</DataArray>" << endl;
}

void ElasticProperties::WriteStressIncrementsVTK(const int tStep) const
{
    stringstream buffer;
    std::vector<int> DataTypes {PDScalars};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        VTK::WriteSymmetricTensor(buffer, StressIncrements, "Sigma");

        buffer << "<DataArray type = \"Float64\" Name = \"" << "Pressure" <<
                    "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < StressIncrements.sizeZ(); ++k)
        for(int j = 0; j < StressIncrements.sizeY(); ++j)
        for(int i = 0; i < StressIncrements.sizeX(); ++i)
        {
            buffer << StressIncrements(i,j,k).Pressure() << endl;
        }
        buffer << "</DataArray>" << endl;

        buffer << "<DataArray type = \"Float64\" Name = \"" << "vMises" <<
                    "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < StressIncrements.sizeZ(); ++k)
        for(int j = 0; j < StressIncrements.sizeY(); ++j)
        for(int i = 0; i < StressIncrements.sizeX(); ++i)
        {
            buffer << Stresses(i,j,k).Mises() << endl;
        }
        buffer << "</DataArray>" << endl;
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "StressIncrements", tStep);
}


void ElasticProperties::WriteVelocityGradientVTK(int tStep)
{

    stringstream buffer;
    std::vector<int> DataTypes {PDScalars};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        VTK::WriteMatrix(buffer, VelocityGradient, "L");
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "VelocityGradient", tStep);
}

void ElasticProperties::PrintPointStatistics(const int x, const int y, const int z) const
{
    Info::WriteCoordinate(x, y, z, dx);
    stringstream outstream;
    outstream << "Stresses: " << Stresses(x, y, z).print() << endl;
    outstream << "Strains: " << Strains(x, y, z).print() << endl;
    outstream << "Effective EigenStrains: " << EffectiveEigenStrains(x, y, z).print() << endl;
    outstream << "Effective ElasticConstants: " << endl
              << EffectiveElasticConstants(x, y, z).print() << endl;

    Info::WriteSimple(outstream.str());
}

vStress ElasticProperties::WriteStressStrainData(std::string filename, std::string LDflag)
{
    vStress aveStress;
    aveStress.set_to_zero();

    STORAGE_LOOP_BEGIN(i,j,k,Stresses,0)
    {
        aveStress += Stresses(i,j,k);
    }
    STORAGE_LOOP_END
    aveStress /= (Nx*Ny*Nz);

    ofstream outputFile;
    outputFile.open(filename, ios::app|ios::out);
    outputFile << AppliedStrain[0]  << ", " << aveStress[0] << ", " <<
                  AppliedStrain[1]  << ", " << aveStress[1] << ", " <<
                  AppliedStrain[2]  << ", " << aveStress[2] << ", " <<
                  AppliedStrain[3]  << ", " << aveStress[3] << ", " <<
                  AppliedStrain[4]  << ", " << aveStress[4] << ", " <<
                  AppliedStrain[5]  << ", " << aveStress[5] << ", " <<
                                               aveStress.Mises() << endl;
    outputFile.close();
    return aveStress;
}

void ElasticProperties::SetVelocityGradient(double* (&rlDefGrad)[9], double dt)
{
    if(VelocityGradient.IsAllocated())
    {
        double velnorm = 1.0/dt;

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,VelocityGradient,0,)
        {
            VelocityGradient(i,j,k)(0,0) = rlDefGrad[0][k + Nz*(j + Ny*i)]*velnorm;
            VelocityGradient(i,j,k)(0,1) = rlDefGrad[1][k + Nz*(j + Ny*i)]*velnorm;
            VelocityGradient(i,j,k)(0,2) = rlDefGrad[2][k + Nz*(j + Ny*i)]*velnorm;
            VelocityGradient(i,j,k)(1,0) = rlDefGrad[3][k + Nz*(j + Ny*i)]*velnorm;
            VelocityGradient(i,j,k)(1,1) = rlDefGrad[4][k + Nz*(j + Ny*i)]*velnorm;
            VelocityGradient(i,j,k)(1,2) = rlDefGrad[5][k + Nz*(j + Ny*i)]*velnorm;
            VelocityGradient(i,j,k)(2,0) = rlDefGrad[6][k + Nz*(j + Ny*i)]*velnorm;
            VelocityGradient(i,j,k)(2,1) = rlDefGrad[7][k + Nz*(j + Ny*i)]*velnorm;
            VelocityGradient(i,j,k)(2,2) = rlDefGrad[8][k + Nz*(j + Ny*i)]*velnorm;
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    else
    {
        Info::WriteExit("VelocityGradient is not allocated!", "ElasticProperties",
                                                       "SetVelocityGradient()");
        exit(13);
    }
}

void ElasticProperties::SetRotations(Orientations& OR, double dt)
{
    if(VelocityGradient.IsAllocated())
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,VelocityGradient,0,)
        {
            // Update rotation via Euler-Rodrigues formula
            dMatrix3x3 W = (VelocityGradient(i,j,k) - VelocityGradient(i,j,k).transposed())*0.5;
//                cout<<EP.VelocityGradient(i,j,k).print()<<endl;
//                cout<<EP.VelocityGradient(i,j,k).transposed().print()<<endl;
//                cout << W.print() << endl;
            dMatrix3x3 RotMat;
            RotMat.set_to_unity();

            double wnorm = 0.0;
            for(int ii = 0; ii < 3; ii++)
            for(int jj = 0; jj < 3; jj++)
            {
                wnorm += W(ii,jj)*W(ii,jj);
            }
            wnorm = sqrt(wnorm*0.5);

            if (wnorm > 0.0)
            {
                RotMat += W*(sin(wnorm*dt)/wnorm)+W*W*((1.0-cos(wnorm*dt))/(wnorm*wnorm)); //Euler-Rodrigues formula
            }

            OR.Quaternions(i,j,k).set(RotMat*OR.Quaternions(i,j,k).RotationMatrix);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    else
    {
        Info::WriteExit("VelocityGradient or OR.Rotations is not allocated!",
                                           "ElasticProperties", "SetRotations()");
        exit(13);
    }
}

void ElasticProperties::SetStressesRX(PhaseField& Phase, PlasticProperties& PFP, BoundaryConditions& BC, std::vector<int> targetPhases)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Stresses,0,)
    {
        if(Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin();
                     alpha != Phase.Fields(i, j, k).cend(); ++alpha)
            {
                int phaseIndex = alpha->index;
                int thPhaseIndex = Phase.FieldsStatistics[phaseIndex].Phase;
                for (auto jt = targetPhases.begin(); jt < targetPhases.end(); jt++)
                {
                    if (*jt == phaseIndex)
                    {
                        vStress temp = Stresses(i,j,k) * alpha->value;
                        Stresses(i,j,k) *= (1 - alpha->value);
                        vStrain tempStrain = EffectiveElasticConstants(i,j,k).inverted()*temp;
                        temp = PhaseElasticConstants[thPhaseIndex] * tempStrain;
//                        double hydroStress = 0.0;
                        for(int kt = 0; kt < 3; kt++)
                        {
//                            hydroStress += temp[kt];
                            temp[kt+3] = 0.0;
                        }
//                        hydroStress /= 3;
//                        for(int kt = 0; kt < 3; kt++)
//                        {
//                            temp[kt] = hydroStress;
//                        }
//                        vStress temp = EffectiveElasticConstants(i,j,k) * (EigenStrains[thPhaseIndex]) * alpha->value;
                        Stresses(i,j,k) += temp * 0.1;
                    }
                }
            } // end alpha
        }
        else
        {
            int locIndex = Phase.Fields(i,j,k).front().index;
            int thPhaseIndex = Phase.FieldsStatistics[locIndex].Phase;
            for (auto jt = targetPhases.begin(); jt < targetPhases.end(); jt++)
            {
                if (*jt == locIndex)
                {
                    vStress temp = Stresses(i,j,k);
//                    Stresses(i,j,k) *= (1 - alpha->value);
                    vStrain tempStrain = EffectiveElasticConstants(i,j,k).inverted()*temp;
                    temp = PhaseElasticConstants[thPhaseIndex] * tempStrain;
//                    double hydroStress = 0.0;
                    for(int kt = 0; kt < 3; kt++)
                    {
//                        hydroStress += temp[kt];
                        temp[kt+3] = 0.0;
                    }
//                    hydroStress /= 3;
//                    for(int kt = 0; kt < 3; kt++)
//                    {
//                        temp[kt] *= 0.5;
////                        temp[kt] = hydroStress;
//                    }
//                    vStress temp = EffectiveElasticConstants(i,j,k) * (EffectiveEigenStrains(i,j,k));
                    Stresses(i,j,k) = temp * 0.1;
                }
            }
        } // end interface
    } // end loop space
    OMP_PARALLEL_STORAGE_LOOP_END
    SetBoundaryConditions(BC);
}

void ElasticProperties::SmoothEigenstrain()
{
    Storage3D<vStrain, 0> EigenStrainsTmp;

//    const double EigenStrainStencil[3][3][3] = {{{1.0/64.0,   1.0/32.0, 1.0/64.0},
//                                                 {1.0/32.0,   1.0/16.0, 1.0/32.0},
//                                                 {1.0/64.0,   1.0/32.0, 1.0/64.0}},
//
//                                                {{1.0/32.0,   1.0/16.0, 1.0/32.0},
//                                                 {1.0/16.0,   1.0/8.0, 1.0/16.0},
//                                                 {1.0/32.0,   1.0/16.0, 1.0/32.0}},
//
//                                                {{1.0/64.0,   1.0/32.0, 1.0/64.0},
//                                                 {1.0/32.0,   1.0/16.0, 1.0/32.0},
//                                                 {1.0/64.0,   1.0/32.0, 1.0/64.0}}};
    double d100 = 1.0/2.0;
    double d110 = 1.0/(3.0 * sqrt(2.0));
    double d111 = 1.0/(6.0 * sqrt(3.0));

    const double EigenStrainStencil[3][3][3] = {{{d111,   d110, d111},
                                                 {d110,   d100, d110},
                                                 {d111,   d110, d111}},

                                                {{d110,   d100, d110},
                                                 {d100,   1.0, d100},
                                                 {d110,   d100, d110}},

                                                {{d111,   d110, d111},
                                                 {d110,   d100, d110},
                                                 {d111,   d110, d111}}};

//    const double EigenStrainStencil[3][3][3] = {{{1.0,   1.0, 1.0},
//                                                 {1.0,   1.0, 1.0},
//                                                 {1.0,   1.0, 1.0}},
//
//                                                {{1.0,   1.0, 1.0},
//                                                 {1.0,   1.0, 1.0},
//                                                 {1.0,   1.0, 1.0}},
//
//                                                {{1.0,   1.0, 1.0},
//                                                 {1.0,   1.0, 1.0},
//                                                 {1.0,   1.0, 1.0}}};

//    const double EigenStrainStencil[3][3][3] = {{{1.0/30.0,   1.0/10.0, 1.0/30.0},
//                                              {1.0/10.0,   7.0/15.0, 1.0/10.0},
//                                              {1.0/30.0,   1.0/10.0, 1.0/30.0}},
//
//                                             {{1.0/10.0,   7.0/15.0, 1.0/10.0},
//                                              {7.0/15.0, -64.0/15.0, 7.0/15.0},
//                                              {1.0/10.0,   7.0/15.0, 1.0/10.0}},
//
//                                             {{1.0/30.0,   1.0/10.0, 1.0/30.0},
//                                              {1.0/10.0,   7.0/15.0, 1.0/10.0},
//                                              {1.0/30.0,   1.0/10.0, 1.0/30.0}}};


//    CreepProperties& CP         = *(static_cast<CreepProperties*>(findOPObject(Objects, "CreepProperties", thisclassname, "Solve", true)));

//    int Nx = CP.Nx;
//    int Ny = CP.Ny;
//    int Nz = CP.Nz;

    EigenStrainsTmp.Allocate(Nx, Ny, Nz, 1);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EigenStrainsTmp,EigenStrainsTmp.Bcells(),)
    {
        EigenStrainsTmp(i,j,k).set_to_zero();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EigenStrainsTmp,0,)
    {
        double WeightSum = 0;
        for(int a = -1; a <= +1; a++)
        for(int b = -1; b <= +1; b++)
        for(int c = -1; c <= +1; c++)
        {
            if (i+a < Nx and j+b < Ny and k+c < Nz and i+a >= 0 and j+b >= 0 and k+c >= 0)
            {
                EigenStrainsTmp(i,j,k) += EffectiveEigenStrains(i+a,j+b,k+c) * EigenStrainStencil[a+1][b+1][c+1];
                WeightSum += EigenStrainStencil[a+1][b+1][c+1];
//              cout << "CreepStrain " << CP.EffectiveCreepStrains(i+a,j+b,k+c).print() << "  Stencil " << CreepStrainStencil[a+1][b+1][c+1] << endl;
            }
        }
//      if (k == 0)
//      cout << i << "," << j << "," << k << "  WeightSum " << WeightSum << "  Eigenstrains " << EigenStrainsTmp(i,j,k).print()<< endl;
        if (fabs(WeightSum) > 0)
        {
//          cout << "CreepStrain " << CreepStrainTmp(i,j,k).print() << endl;
            EigenStrainsTmp(i,j,k) /= WeightSum;
//          cout << "CreepStrain " << CreepStrainTmp(i,j,k).print() << "  WeightSum " << WeightSum << endl;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    EffectiveEigenStrains = EigenStrainsTmp;
}

void ElasticProperties::ApplyDamage(DamageProperties& DP)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EffectiveElasticConstants, EffectiveElasticConstants.Bcells(), )
    {
        EffectiveElasticConstants(i, j, k) *= (1.0 - DP.EffectiveDamage(i, j, k));
#ifdef OMP
#pragma omp critical
#endif
        if (i >= 0 and i < Nx and j >= 0 and j < Ny and k >= 0 and k < Nz)
        {
            AverageElasticConstants += EffectiveElasticConstants(i, j, k);
            AverageCompliences += EffectiveElasticConstants(i, j, k).inverted();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    AverageElasticConstants /= Nx*Ny*Nz;
    AverageCompliences /= Nx*Ny*Nz;
}

void ElasticProperties::ApplyRotation(Orientations& OR)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, EffectiveElasticConstants, EffectiveElasticConstants.Bcells(), )
    {
		if (abs(OR.Quaternions(i, j, k).RotationMatrix.trace() - 3.0) > 1.0e-4)
		{
			EffectiveElasticConstants(i, j, k).rotate(OR.Quaternions(i, j, k).RotationMatrix);
			EffectiveEigenStrains(i, j, k).rotate(OR.Quaternions(i, j, k).RotationMatrix);
		}
#ifdef OMP
#pragma omp critical
#endif
        if (i >= 0 and i < Nx and j >= 0 and j < Ny and k >= 0 and k < Nz)
        {
            AverageElasticConstants += EffectiveElasticConstants(i, j, k);
            AverageCompliences += EffectiveElasticConstants(i, j, k).inverted();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    AverageElasticConstants /= Nx*Ny*Nz;
    AverageCompliences /= Nx*Ny*Nz;
}


double ElasticProperties::Energy(void)
{
    double Energy = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Strains, 0, reduction(+:Energy))
    {
        Energy += PointEnergy(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    return Energy*dx*dx*dx;
}

double ElasticProperties::AverageEnergyDensity(void)
{
    double Energy = 0.0;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Strains, 0, reduction(+:Energy))
    {
        Energy += PointEnergy(i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    return Energy/double(Nx*Ny*Nz);
}

void ElasticProperties::CalculateInterfaceEnergyContribution(PhaseField& Phase, InterfaceEnergy& IE)
{
	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Strains, 0,)
	{
		if(Phase.Fields(i,j,k).flag)
		{
			for(auto alpha  = Phase.Fields(i, j, k).cbegin();
					 alpha != Phase.Fields(i, j, k).cend() - 1; ++alpha)
			for(auto  beta  = alpha + 1;
					  beta != Phase.Fields(i, j, k).cend(); ++beta)
			{
				IE.RawIntEnergy3D(i,j,k).add_sym(alpha->index, beta->index, PointEnergy(i,j,k)*Phase.Eta*10.0);
			}
		}
	}
	OMP_PARALLEL_STORAGE_LOOP_END
}

double ElasticProperties::PointEnergy(int i, int j, int k)
{
    vStrain elStrain;
    elStrain.set_to_zero();
    dMatrix6x6 EffectiveCompliance = EffectiveElasticConstants(i,j,k).inverted();

    for(int ii = 0; ii < 6; ii++)
    for(int jj = 0; jj < 6; jj++)
    {
        elStrain[ii] += Stresses(i,j,k)[jj]*EffectiveCompliance(ii,jj);
    }

    double Energy = 0.0;

    for(int ii = 0; ii < 6; ii++)
    for(int jj = 0; jj < 6; jj++)
    {
        Energy += elStrain[ii]*EffectiveElasticConstants(i,j,k)(ii,jj)*elStrain[jj];
    }
    return Energy*0.5;
}

void ElasticProperties::WriteEnergyVTK(int tStep)
{
    stringstream buffer;
    std::vector<int> DataTypes {PDScalars};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        buffer << "<DataArray type = \"Float64\" Name = \"" << "ElasticEnergy" <<
                    "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            buffer << PointEnergy(i,j,k) << endl;
        }
        buffer << "</DataArray>" << endl;
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "ElasticEnergy", tStep);
}  //  WriteVTK
}// namespace openphase
