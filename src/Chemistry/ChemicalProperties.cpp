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

#include "Info.h"
#include "PhaseField.h"
#include "Tools/Node.h"
#include "BoundaryConditions.h"
#include "Tools/UserInterface.h"
#include "VTK.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "Chemistry/ChemicalProperties.h"
#include "Chemistry/PeriodicTable.h"
#include "Compositions.h"

//TODO: use internal/external loops instead of for()for()for()
//TODO: full documentation of the functions in /***/
//TODO: use static functions which are initialised instead of functions that
//      have to be called all the time int Nsites <-> int Nsites()
//TODO: VTK/RawData is ok?

namespace opensim
{

using namespace std;

void ChemicalProperties::Initialize(Settings& locSettings,
                                    const int boundarysize)
{
    /**This function will take system properties that have been read in before,
    simulation box dimensions, number of components, number of thermodynamic
    phases and will initialize the storages. This will also initialize the
    phase-storage, which by itself will initialize the sublattice-storage.
    (ReadManualInput() or ReadInput()+TQ.ReadInput() has to be called
    before!).*/

    thisclassname = "ChemicalProperties";
    Nx      = locSettings.Nx;
    Ny      = locSettings.Ny;
    Nz      = locSettings.Nz;
    AtStart = true;
    Threshold = DBL_EPSILON;
    Eta        = locSettings.Eta;
    dx = locSettings.dx;
    dt = locSettings.dt;

    Ncomp = getNcomp();
    Nphases = getNphases();

    PhaseInput phinput;
    phinput.Components = Component;
    phinput.RefElements = RefElements;

    for(int i = 0; i < Nphases; i++)
    {
        phinput.Number = i;
        Phase[i].Initialize(locSettings, boundarysize, phinput);
    }

    Ncons = getNcons();
    Nsubs = getNsubs();
    TotalNsubs = getTotalNsubs();
    Nsites = getNsites();
    SublIdx = getSublIdx();
    ConsIdx = getConsIdx();
    NsitesWithI = getNsitesWithI();
    SortedElementMatrix = CalculateSortedElementMatrix();

    TotInitial.resize(Ncomp);

    CalculateMolefractionLimits();

    if(locSettings.CP.diffMod == ChemicalProperties::DiffusionModel::EquilibriumPartitioning)
    {
        for(int i = 0; i < Nphases; i++)
        {
            Phase[i].StoichiometricFlag = false;
        }
    }
    else if(locSettings.CP.diffMod == ChemicalProperties::DiffusionModel::FiniteInterfaceDissipation)
    {
        for(int i = 0; i < Nphases; i++)
        {
            Phase[i].StoichiometricFlag = Phase[i].get_isStoichiometric();
        }
    }

    PrintData();

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void ChemicalProperties::InitializeTracer(const int boundarysize)
{
    TotalTracer.Allocate(Nx, Ny, Nz, {NumIsotopes}, boundarysize);
    PhaseDotTracer.Allocate(Nx, Ny, Nz, {NumIsotopes,3,3,3}, boundarysize);
    minitialTracer.resize(NumIsotopes);
}

void ChemicalProperties::ReadInput(Settings& locSettings, string FileName)
{
    /**Read number of components, their names and atomic weights, number of
    phases and also their names. This function will not read the sublattice
    model also present in the same input file. This will be the base data
    needed to read the TDB/GES5/manual input files.*/

    thisclassname = "ChemicalProperties";
    PeriodicTable PT;

    fstream inp(FileName.c_str(), ios::in);

    if (!inp)
    {
        string thisclassname = "ChemicalProperties";
        Info::WriteExit("File \"" + FileName + "\" could not be opened",
                                                thisclassname, "ReadInput");
        exit(1);
    };

    Info::WriteLine();
    Info::WriteLineInsert("ChemicalProperties");
    Info::WriteStandard("Source", FileName);

    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

	string tmpDiff = UserInterface::ReadParameterS(inp, moduleLocation, "DiffusionModel", false, "NONE");
	std::transform(tmpDiff.begin(), tmpDiff.end(), tmpDiff.begin(), ::toupper);
	tmpDiff.erase(std::remove(tmpDiff.begin(), tmpDiff.end(), ' '), tmpDiff.end());
	if (tmpDiff == "NONE")
	{
		diffMod = DiffusionModel::None;
	}
	if (tmpDiff == "EQUILIBRIUMPARTITIONING" or tmpDiff == "EQUILIBRIUMPARTITION" or tmpDiff == "EQUILIBRIUMPARTITIONDIFFUSION")
	{
		diffMod = DiffusionModel::EquilibriumPartitioning;
	}
	if (tmpDiff == "SUBLATTICEMODEL" or tmpDiff == "SUBLATTICE" or tmpDiff == "SUBLATTICEDIFFUSION" or tmpDiff == "FINITEINTERFACEDISSIPATION")
	{
		diffMod = DiffusionModel::FiniteInterfaceDissipation;
	}

    bool endofnames = false;
    int n = 0;

    while(!(endofnames))
    {
        stringstream converter;
        converter << string("Comp_") << n;
        if(UserInterface::FindParameter(inp, moduleLocation, converter.str()) != -1)
        {
            string tmp = UserInterface::ReadParameterF(inp, moduleLocation, converter.str());
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
            if(tmp == "VA")
            {
                cerr << "\"VA\" as a component is not allowed in "
                     << "ChemicalProperties.opi" << endl;
                //exit(0);
            }
            Elementnames.push_back(tmp);
            n++;
        }
        else
        {
            endofnames = true;
        }
    }

    for(unsigned int i = 0; i < Elementnames.size(); i++)
    {
        Element tempComp;
        tempComp.Index = i;
        tempComp.Value = 0.0;
        tempComp.Min = 0.0;
        tempComp.Max = 1.0;
        tempComp.Name = Elementnames[i];
        tempComp.isVacancy = false;
        tempComp.isStoichiometric = false;
        tempComp.Major = false;
        Component.push_back(tempComp);
    }

    endofnames = false;
    n = 0;

    while(!(endofnames))
    {
        stringstream converter;
        converter << string("Phase_") << n;
        if(UserInterface::FindParameter(inp, moduleLocation, converter.str()) != -1)
        {
            string tmp = UserInterface::ReadParameterF(inp, moduleLocation, converter.str());
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
            Phasenames.push_back(tmp);
            n++;
        }
        else
        {
            endofnames = true;
        }
    }

    for(unsigned int n = 0; n < Phasenames.size(); n++)
    {
        ThermodynamicPhase tmpPhase;
        tmpPhase.Name = Phasenames[n];
        //tmpPhase.Index = n;
        Phase.push_back(tmpPhase);
    }

    if(Component.size() != 0)
    for (unsigned int n = 0; n < Phasenames.size(); n++)
    {
        stringstream converter;
        converter << string("RefElement_") << n;
        if (UserInterface::FindParameter(inp, moduleLocation, converter.str()))
        {
            string tmp = UserInterface::ReadParameterF(inp, moduleLocation, converter.str(),false,Elementnames[0]);
            std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::toupper);
            RefElements.push_back(tmp);
        }
    }

    int moduleLocationTracer =UserInterface::FindModuleLocation(inp, "TracerDiffusion");
    if(moduleLocationTracer != 0 && UserInterface::FindParameter(inp, moduleLocationTracer, string("Isotope")) != -1)
    {
    	Isotope = UserInterface::ReadParameterB(inp, moduleLocationTracer, string("Isotope"));
    	ReadInputTracer(inp, moduleLocationTracer);
    }
    else
    {
    	Isotope = false;
    }

    inp.close();
    Info::WriteLine();

	Initialize(locSettings, 0);
}


void ChemicalProperties::ReadInputTracer(fstream& inp, int moduleLocation)
{
	int Ncomp = Component.size();
    for(int comp = 0; comp < Ncomp; comp++)
    {
        stringstream converter;
        converter << string("C0_") << Component[comp].Name << string("_Tracer_Amount");
        if(UserInterface::FindParameter(inp, moduleLocation, converter.str()) != -1)
        {
            IsotopeDefined.push_back(Component[comp].Index);

        }
    }
    NumIsotopes = IsotopeDefined.size();
    InitialTracerAmount.resize(NumIsotopes);
    InitialTracerPosition.resize(NumIsotopes);

    //ReadIsotopeConcentration
    for(int i = 0; i <NumIsotopes; i++)
    {
        stringstream converter;
        converter << string("C0_") << Component[IsotopeDefined[i]].Name << string("_Tracer_Amount");
        InitialTracerAmount[i]= UserInterface::ReadParameterD(inp, moduleLocation, converter.str());
    }

    //ReadIsotopePosition
    for(int i = 0; i <NumIsotopes; i++)
    {
        stringstream converter;
        converter << string("C0_") << Component[IsotopeDefined[i]].Name << string("_Tracer_X");
        InitialTracerPosition[i]= UserInterface::ReadParameterD(inp, moduleLocation, converter.str());
    }
}

vector<int> ChemicalProperties::getNcons(void)
{
    /**This function will return a vector with the total amount of constituents
    for each phase.*/
    std::vector<int> tempNcons(Nphases);
    for(int i = 0; i < Nphases; i++)
    {
        tempNcons[i] = Phase[i].Ncons;
    }
    return tempNcons;
}

void ChemicalProperties::ReadCompositionInput(string FileName)
{
    /**This function will read the composition input file and write this data
    to a vector in ThermodynamicPhase.*/

    fstream inp(FileName.c_str(), ios::in);

    if (!inp)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be opened",
                                  thisclassname, "ReadCompositionInput");
        exit(1);
    };

    Info::WriteLine();
    Info::WriteLineInsert("Composition properties");
    Info::WriteStandard("Source", FileName);

    int moduleLocation = UserInterface::FindModuleLocation(inp, "Composition");

    /* Composition Read Input */

    for(int alpha = 0; alpha < Nphases; alpha++)
    {
        int cons = 0;
        for(int s = 0; s < Phase[alpha].Nsubs; s++)
        {
            double unity = 0;
            for(int n = 0; n < Phase[alpha].Sublattice[s].Ncons; n++)
            {
                stringstream converter;
                converter << string("C0_") << Phase[alpha].Name << string("_")
                          << s << string("_") <<
                                 Phase[alpha].Sublattice[s].Constituent[n].Name;

                Phase[alpha].Initial({cons}) =
                            UserInterface::ReadParameterD(inp, moduleLocation, converter.str());
                unity += Phase[alpha].Initial({cons});
                cons++;
            }
            if(unity < 1.0-1E-6 or unity > 1.0+1E-6)
            {
                cerr << "Site fractions don't add up to unity, Phase "
                     << Phase[alpha].Name << ", Sublattice " << s
                     << ". Unity = " << unity << endl;
                //exit(1);
            }
        }
    }

    /* Second Composition Read Input (for composition slopes)*/

    for(int alpha = 0; alpha < Nphases; alpha++)
    {
        int cons = 0;
        for(int s = 0; s < Phase[alpha].Nsubs; s++)
        {
            double unity = 0;
            for(int n = 0; n < Phase[alpha].Sublattice[s].Ncons; n++)
            {
                stringstream converter;
                converter << string("C1_") << Phase[alpha].Name << string("_")
                          << s << string("_") <<
                                 Phase[alpha].Sublattice[s].Constituent[n].Name;

                Phase[alpha].Initial2({cons})
                = UserInterface::ReadParameterD(inp, moduleLocation,
                                                converter.str(), false,
                                                Phase[alpha].Initial({cons}));

                unity += Phase[alpha].Initial2({cons});
                cons++;
            }
            if(unity < 1.0-1E-6 or unity > 1.0+1E-6)
            {
                cerr << "Site fractions don't add up to unity, Phase "
                     << Phase[alpha].Name << ", Sublattice " << s
                     << ". Unity = " << unity << endl;
                //exit(1);
            }
        }
    }

    inp.close();
    Info::WriteLine();
}

void ChemicalProperties::PrintData(void)
{
    /**This function will print the names of all components stored in the local
    copy of the thermodynamic system, the index of each component as stored in
    the TDB/GES5/manual input file, the names of all thermodynamic phases and
    its index, as well as the sublattice model. Vacancies are represented with
    a -1.*/

    Info::WriteLineInsert("Data stored in " + thisclassname +
                                               ". Position in TDB in brackets");

    string ctmp;

    for(int i = 0; i < Ncomp; i++)
    {
        //ctmp += "[" + to_string(Component[i].Index) + "] " + Component[i].Name;
        ctmp += "[] " + Component[i].Name;
        if(i < Ncomp-1)
        ctmp += ", ";
    }

    Info::WriteStandard("Components in the system", ctmp);

    string ptmp;

    for(int i = 0; i < Nphases; i++)
    {
        //ptmp += "[" + to_string(Phase[i].Index) + "] " + Phase[i].Name;
        ptmp += "[] " + Phase[i].Name;

        string tmp;

        for(int j = 0; j < Phase[i].Nsubs; j++)
        {
            std::ostringstream ostrs;
            tmp.append(" (");
            for(int k = 0; k < Phase[i].Sublattice[j].Ncons; k++)
            {
                tmp.append(Phase[i].Sublattice[j].Constituent[k].Name);
                if(k < Phase[i].Sublattice[j].Ncons-1)
                tmp.append(",");
            }
            tmp.append(")_");
            ostrs << Phase[i].Sublattice[j].Site;
            tmp.append(ostrs.str());
        }

        ptmp += tmp;

        if(i < Nphases-1)
        ptmp += ", ";
    }

    Info::WriteStandard("Phases in the system", ptmp);
    Info::WriteLine();
}

void ChemicalProperties::CalculateMoleFractions(Composition& Cx)
{
    /**This function will cal*/
    for(int n = 0; n < Nphases; n++)
    {
        Phase[n].CalculateMoleFractions(Cx);
    }
}

vector<double> ChemicalProperties::getTotalComposition(Composition& Cx, PhaseField& Phi, int x, int y, int z)
{
    //TODO please use Composition::getTotalComposition
    return Cx.getTotalComposition(Phi,x,y,z);
}

vector<double> ChemicalProperties::getWeightPercent(Composition& Cx, PhaseField& Phi, int x, int y, int z)
{
    //TODO please use Composition::getWeightPercent

    /**This function will loop over all components and calculate
    the local weight percent of each element, by multiplying the total mole
    fractions with the atomic mass.*/

    PeriodicTable PT;
    vector<double> result(Ncomp, 0.0);
    vector<double> AtomicWeight(Ncomp, 0.0);
    vector<double> Total = getTotalComposition(Cx,Phi,x,y,z);

    for(int comp = 0; comp < Ncomp; comp++)
    AtomicWeight[comp] = PT.GetData(Component[comp].Name).AtomicWeight;

    double div = 0.0;
    for(int comp = 0; comp < Ncomp; comp++)
    {
        div += Total[comp]*AtomicWeight[comp];
    }
    for(int comp = 0; comp < Ncomp; comp++)
    {
        result[comp] = Total[comp]*AtomicWeight[comp]/div;
    }

    return result;
}

void ChemicalProperties::Remesh(int newNx, int newNy, int newNz,
                                                         BoundaryConditions& BC)
{
    /**This function will change the dimension of the simulation box. Nx, Ny
    and Nz will be overwritten, the 3D storage Total will be remeshed as well as
    the dimension of the Composition and the Mole fraction storage of each
    phase.*/

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    for(int i = 0; i < Nphases; i++)
    {
        Phase[i].Nx = newNx;
        Phase[i].Ny = newNy;
        Phase[i].Nz = newNz;
        //Phase[i].Composition.Remesh(newNx, newNy, newNz);
        //Phase[i].MoleFractions.Remesh(newNx, newNy, newNz);
    }
    SetBoundaryConditions(BC);
    Info::WriteStandard(thisclassname, "Remeshed");
}

void ChemicalProperties::MoveFrame(int dx, int dy, int dz, 
                                                         BoundaryConditions& BC)
{
    /**This function will move the 3D storages Total, Composition and
    MoleFractions, when given an increment in each direction.*/

    int xBeg = (dx >= 0) + (dx < 0)*(Nx) - 1;
    int xEnd = (dx >= 0)*(Nx) + (dx < 0) - 1;
    int xInc = 1 - 2*(dx < 0);

    int yBeg = (dy >= 0) + (dy < 0)*(Ny) - 1;
    int yEnd = (dy >= 0)*(Ny) + (dy < 0) - 1;
    int yInc = 1 - 2*(dy < 0);

    int zBeg = (dz >= 0) + (dz < 0)*(Nz) - 1;
    int zEnd = (dz >= 0)*(Nz) + (dz < 0) - 1;
    int zInc = 1 - 2*(dz < 0);

    SetBoundaryConditions(BC);

    for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
    {
        for(int n = 0; n < Nphases; n++)
        {
            //Phase[n].Composition(i,j,k) =
            //                       Phase[n].Composition(i + dx, j + dy, k + dz);
            //Phase[n].MoleFractions(i,j,k) =
            //                     Phase[n].MoleFractions(i + dx, j + dy, k + dz);
        }
    }

    SetBoundaryConditions(BC);
    Info::WriteStandard(thisclassname, "Frame moved");
}

void ChemicalProperties::SetBoundaryConditions(BoundaryConditions& BC)
{
    /**Sets boundary conditions for the different composition storages of every
    thermodynamic phase.*/
    //TODO: Boundary conditions for Total needed? Yes :)
    /*BC.SetX(TotalTracer);
    BC.SetY(TotalTracer);
    BC.SetZ(TotalTracer);*/
}

void ChemicalProperties::SetBoundaryConditionsTracer(BoundaryConditions& BC)
{
    /**Sets boundary conditions for the tracer composition.*/
    BC.SetX(TotalTracer);
    BC.SetY(TotalTracer);
    BC.SetZ(TotalTracer);
}

int ChemicalProperties::PhaseNumber(string name)
{
    /**This function is used to address thermodynamic phases by their names. It
    returns the index, given the phase name as a string.*/
    transform(name.begin(), name.end(), name.begin(), ::toupper);
    int phasenr = -1;
    for(int i = 0; i < Nphases; i++)
    if(name == Phase[i].Name)
    {
        phasenr = i;
    }
    if(phasenr < 0)
    {
        cerr << "Error in ChemicalProperties::PhaseNumber(), no phase available"
             << " with the name of " << name << endl;
        //exit(1);
    }
    return phasenr;
}

int ChemicalProperties::getNphases(void)
{
    /**This function will  the number of stored thermodynamic phases in the
    system.*/
    return Phase.size();    
}
int ChemicalProperties::getNcomp(void)
{
    /**This function will return the number of stored components in the
    system.*/
    return Component.size();
}

vector<double> ChemicalProperties::getNsites(void)
{
    /**This function will loop over all phases and return the total number of
    sites for each phase.*/
    int locNphases = Nphases;
    vector<double> temp(locNphases, 0.0);
    for(int n = 0; n < locNphases; n++)
    {
        temp[n] = Phase[n].Nsites;
    }
    return temp;
}

vector<int> ChemicalProperties::getNsubs(void)
{
    /**This function will return a vector with the number of sublattices for
    each thermodynamic phase.*/
    int locNphases = Nphases;
    vector<int> temp(locNphases, 0);
    for(int n = 0; n < locNphases; n++)
    {
        temp[n] = Phase[n].Nsubs;
    }
    return temp;
}

int ChemicalProperties::getTotalNsubs(void)
{
    /**This function will return an integer with the number of all sublattices
    in the thermodynamic system.*/
    int locNphases = Nphases;
    int temp = 0;
    for(int n = 0; n < locNphases; n++)
    {
        temp += Phase[n].Nsubs;
    }
    return temp;
}

vector<vector<int> > ChemicalProperties::getSublIdx(void)
{
    /**This function will return a vector of a vector, containing the first and
    the last constituent by its index, to loop over all constituents in each
    sublattice. For more information, see ThermodynamicPhase:SublIdx().*/
    int locNphases = Nphases;
    vector<vector<int> > temp;
//    vector<vector<int> > temp(locNphases, {0});
    for(int n = 0; n < locNphases; n++)
    {
    	temp.push_back(Phase[n].SublIdx);
//        temp[n] = Phase[n].SublIdx;
    }
    return temp;
}

vector<vector<int> > ChemicalProperties::getConsIdx(void)
{
    /**This function will return a vector of a vector, containing the index of
    each constituent in its order given by the sublattice model.*/
    int locNphases = Nphases;
    vector<vector<int> > temp;
//    vector<vector<int> > temp(locNphases, {0});
    for(int n = 0; n < locNphases; n++)
    {
    	temp.push_back(Phase[n].ConsIdx);
//        temp[n] = Phase[n].ConsIdx;
    }
    return temp;
}

vector<vector<double> > ChemicalProperties::getNsitesWithI(void)
{
    /**This function will return a vector of the output of
    ThermodynamicPhase::NsitesWithI() for each phase. For more information, see
    its description in the ThermodynamicPhase class.*/
    int locNphases = Nphases;
    vector<vector<double> > temp;
//    vector<vector<double> > temp(locNphases, {0.0});
    for(int n = 0; n < locNphases; n++)
    {
    	temp.push_back(Phase[n].NsitesWithI);
//        temp[n] = Phase[n].NsitesWithI;
    }
    return temp;
}

double ChemicalProperties::Nmoles(Composition& Cx, PhaseField& Phi, int x, int y, int z)
{
    /**This function will calculate the total number of moles per point. It will
    loop over all phases, calculate the amount of moles in each phase from
    the stored site fractions and return the sum of all, weighted with the
    phase fractions.*/
    double temp = 0.0;
    for(int alpha = 0; alpha < Nphases; alpha++)
    if(Phi.Fractions(x,y,z)({alpha}) > 0.0)
    {
        temp += Phi.Fractions(x,y,z)({alpha})*Phase[alpha].Nmoles(Cx,x,y,z);
    }
    return temp;
}

void ChemicalProperties::Write(int tStep, bool legacy_format)
{
    /**This function will write the content of */

    /*
    string FileName = UserInterface::MakeFileName(RawDataDir,"Composition_",
                                                              tStep, ".dat");

    ofstream out(FileName.c_str(), ios::out | ios::binary);

    if (!out)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be created",
                                                thisclassname, "Write()");
        exit(1);
    }

    for(int n = 0; n < Nphases; n++)
    {
        STORAGE_LOOP_BEGIN(i,j,k,Phase[n].Composition,0)
        {
            out.write(reinterpret_cast<char*>(Phase[n].Composition(i,j,k).data()),
                             Phase[n].Composition(i,j,k).size()*sizeof(double));
        }
        STORAGE_LOOP_END
    }
    */
}

void ChemicalProperties::Read(BoundaryConditions& BC, PhaseField& Phi,
                              int tStep, bool legacy_format)
{
    /**Read composition in ASCII format stored in the RawData folder. This is
     * used when restarting from a previous calculated time step.*/

    /*
    string FileName = UserInterface::MakeFileName(RawDataDir,"Composition_",
                                                              tStep, ".dat");

    fstream inp(FileName.c_str(), ios::in | ios::binary);

    if (!inp)
    {
        Info::WriteExit("File \"" + FileName + "\" could not be opened",
                                                thisclassname, "Read()");
        exit(1);
    }

    for(int n = 0; n < Nphases; n++)
    {
        STORAGE_LOOP_BEGIN(i,j,k,Phase[n].Composition,0)
        {
            inp.read(reinterpret_cast<char*>(Phase[n].Composition(i,j,k).data()),
                             Phase[n].Composition(i,j,k).size()*sizeof(double));
        }
        STORAGE_LOOP_END
    }

    CalculateMoleFractions();

    for(int x = 0; x < Nx; x++)
    for(int y = 0; y < Ny; y++)
    for(int z = 0; z < Nz; z++)
    {
        vector<double> Total = getTotalComposition(Phi,x,y,z);

        for(int n = 0; n < Ncomp; n++)
        {
            TotInitial[n] += Total[n];
        }
    }
    
    for(int n = 0; n < Ncomp; n++)
    {
        TotInitial[n] /= double(Nx*Ny*Nz);
    }

    SetBoundaryConditions(BC);
    Info::WriteStandard(thisclassname, "Binary input loaded");
    */
}

void ChemicalProperties::WriteVTK(PhaseField& Phi, int tStep)
{
    /**This function will write VTK files for visualization of the 3D storages.
     * Total mole fractions, total weight percent, phase mole fractions as well
     * as site fractions.*/


/*    CalculateMoleFractions();
    stringstream outbufer;

    outbufer << "# vtk DataFile Version 3.0\n";
    outbufer << "Composition\n";
    outbufer << "ASCII\n";
    outbufer << "DATASET STRUCTURED_GRID\n";
    outbufer << "DIMENSIONS " << Nz << " " << Ny << " " << Nx << "\n";
    outbufer << "POINTS " <<  Nx*Ny*Nz << " double\n";
    
    for(int x = 0; x < Nx; x++)
    for(int y = 0; y < Ny; y++)
    for(int z = 0; z < Nz; z++)
    {
        outbufer << x << " " << y << " " << z << "\n";
    }
    
    outbufer << "\n";
    outbufer << "POINT_DATA " << Nx*Ny*Nz << "\n";

    for(int comp = 0; comp < Ncomp; comp++)
    {
        outbufer << "SCALARS TotalCompositionMF_" << Component[comp].Name
                 << " double 1\n";
        outbufer << "LOOKUP_TABLE default\n";

        for(int x = 0; x < Nx; x++)
        for(int y = 0; y < Ny; y++)
        for(int z = 0; z < Nz; z++)
        {
            vector<double> Total = getTotalComposition(Phi,x,y,z);
            outbufer << Total[comp] << "\n";
        }
    }

    for(int comp = 0; comp < Ncomp; comp++)
    {
        outbufer << "SCALARS TotalCompositionWP_" << Component[comp].Name
                 << " double 1\n";
        outbufer << "LOOKUP_TABLE default\n";

        for(int x = 0; x < Nx; x++)
        for(int y = 0; y < Ny; y++)
        for(int z = 0; z < Nz; z++)
        {
            vector<double> WeightPercent(Ncomp, 0.0);
            WeightPercent = getWeightPercent(Phi,x,y,z);

            outbufer << WeightPercent[comp] << "\n";
        }
    }

    for(int alpha = 0; alpha < Nphases; ++alpha)
    for(int comp = 0; comp < Ncomp; comp++)
    {
        outbufer << "SCALARS PhaseCompositionMF_" << Component[comp].Name
                 << "(" << Phase[alpha].Name << ") double 1\n";
        outbufer << "LOOKUP_TABLE default\n";

        STORAGE_LOOP_BEGIN(i,j,k,Phase[alpha].MoleFractions,0)
        {
            if(Phi.Fractions(i,j,k)({alpha}) > 0.0)
            {
                outbufer << Phase[alpha].MoleFractions(i,j,k)({comp}) << "\n";
            }
            else
            {
                outbufer << 0.0 << "\n";
            }
        }
        STORAGE_LOOP_END
    }

    for(int alpha = 0; alpha < Nphases; ++alpha)
    {
        int index = 0;
        for(int sub = 0; sub < Phase[alpha].Nsubs; ++sub)
        for(int cons = 0; cons < Phase[alpha].Sublattice[sub].Ncons; cons++)
        {
            outbufer << "SCALARS PhaseCompositionSF_"
                     << Phase[alpha].Sublattice[sub].Constituent[cons].Name
                     << "#" << sub << "(" << Phase[alpha].Name
                     << ") double 1\n";
            outbufer << "LOOKUP_TABLE default\n";

            STORAGE_LOOP_BEGIN(i,j,k,Phase[alpha].Composition,0)
            {
                if(Phi.Fractions(i,j,k)({alpha}) > 0.0)
                {
                    outbufer << Phase[alpha].Composition(i,j,k)({index}) << "\n";
                }
                else
                {
                    outbufer << 0.0 << "\n";
                }
            }
            STORAGE_LOOP_END
            
            index++;
        }
    }
    string FileName = UserInterface::MakeFileName(VTKDir, "Composition_",
                                                            tStep,".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbufer.rdbuf();
    vtk_file.close();*/

}  //  WriteVTK

void ChemicalProperties::WriteStatistics(Composition& Cx, PhaseField& Phi, int tStep, double dt)
{
    /**obsolete*/

    cout << "WriteStatistics() obsolete, will exit now!" << endl;
    //exit(1);

    vector<double> total(Ncomp, 0);
    vector<double> deviation(Ncomp, 0);
    double sim_time = tStep*dt;

    for(int comp = 0; comp < Ncomp; comp++)
    {
        for(int x = 0; x < Nx; x++)
        for(int y = 0; y < Ny; y++)
        for(int z = 0; z < Nz; z++)
        {
            vector<double> Total = getTotalComposition(Cx,Phi,x,y,z);
            total[comp] += Total[comp];
        }
        
        total[comp] /= double(Nx*Ny*Nz);
        if (AtStart)
        {
            TotInitial[comp] = total[comp];
        }
        deviation[comp] = TotInitial[comp] - total[comp];
    }

    AtStart = false;

    ofstream output_file;
    if (tStep == 0)
    {
        output_file.open("CompositionStatistics.txt", ios::out);
        output_file << "tStep" << "\t\t" << "sim_time" << "\t\t";
        for(int comp = 0; comp < Ncomp; comp++)
        {
            output_file << "total_" << Component[comp].Name << "\t\t"
                        << "deviation_" << Component[comp].Name << "\t\t";
        }
        output_file<< endl;
        output_file.close();
    }

    output_file.open("CompositionStatistics.txt", ios::app);
    output_file << tStep << "\t\t" << sim_time << "\t\t";
    for(int comp = 0; comp < Ncomp; comp++)
    {
        output_file << total[comp]  << "\t\t" << deviation[comp] << "\t\t";
    }
    output_file << endl;
    output_file.close();
}

void ChemicalProperties::CalculateMolefractionLimits(void)
{
    /**This function will analyze the sublattice models of each phase and will
    calculate the minimum and maximum limits for the mole fraction of each
    component.*/
    int locNphases = Nphases;
    int locNcomps  = Ncomp;
    for(int alpha = 0; alpha < locNphases; alpha++)
    for(int comp = 0; comp < locNcomps; comp++)
    {
        Phase[alpha].Component[comp].Min = 0.0;
        Phase[alpha].Component[comp].Max = 1.0;

        int Nsubs = Phase[alpha].Nsubs;
        double totalSitesMAX = 0.0;
        double totalSitesMIN = 0.0;
        double sitesWithCompMAX = 0.0;
        double sitesWithCompMIN = 0.0;
        for(int sub = 0; sub < Nsubs; sub++)
        {
            //COUNT MINIMUM
            if(Phase[alpha].Sublattice[sub].isElementPresent(Component[comp].Index)
            and(Phase[alpha].Sublattice[sub].Ncons < 2))
            {
                sitesWithCompMIN += Phase[alpha].Sublattice[sub].Site;
            }
            if(!(Phase[alpha].Sublattice[sub].hasVacancies
            and(Phase[alpha].Sublattice[sub].Ncons < 2)))
            {
                totalSitesMIN += Phase[alpha].Sublattice[sub].Site;
            }
            //COUNT MAXIMUM
            if(Phase[alpha].Sublattice[sub].isElementPresent(Component[comp].Index))
            {
                sitesWithCompMAX += Phase[alpha].Sublattice[sub].Site;
            }
            if(!(Phase[alpha].Sublattice[sub].hasVacancies)
            and!(Phase[alpha].isInterstitial(comp)))
            {
                totalSitesMAX += Phase[alpha].Sublattice[sub].Site;
            }
            else if(Phase[alpha].isInterstitial(comp))
            {
                totalSitesMAX += Phase[alpha].Sublattice[sub].Site;
            }
        }
        if(totalSitesMIN > 0)
        Phase[alpha].Component[comp].Min = sitesWithCompMIN/totalSitesMIN;
        if(totalSitesMAX > 0)
        Phase[alpha].Component[comp].Max = sitesWithCompMAX/totalSitesMAX;
    }
}
void ChemicalProperties::CalculateTotalMolarVolume(Composition& Cx, PhaseField& Phi)
{
    //TODO: use Composition::CalculateTotalMolarVolume
    CalculateTotalAverage(Cx,Phi);
    TotalMolarVolume = 0.0;
    PeriodicTable PT;

    for(int i = 0; i < Ncomp; i++)
    {
        TotalMolarVolume += PT.GetData(Component[i].Name).MolarVolume
                            *TotalAverage[i];
    }
}
void ChemicalProperties::CalculateTotalAverage(Composition& Cx, PhaseField& Phi)
{
    TotalAverage.resize(Ncomp,0.0);
    int counter = 0;

    for(int x = 0; x < Nx; x++)
    for(int y = 0; y < Ny; y++)
    for(int z = 0; z < Nz; z++)
    {
        counter++;
        vector<double> Total = Cx.getTotalComposition(Phi,x,y,z);

        for(int comp = 0; comp < Ncomp; comp++)
        {
            TotalAverage[comp] += Total[comp];
        }
    }
    
    for(int comp = 0; comp < Ncomp; comp++)
    {
        TotalAverage[comp] /= counter;
    }
}

vector<int> ChemicalProperties::CalculateSortedElementMatrix(void)
{
    vector<int> result(Ncomp);
    vector<string> Names(Ncomp);
    vector<string> SortedNames;

    for(int n = 0; n < Ncomp; n++)
    Names[n] = Component[n].Name;

    SortedNames = Names;
    sort(SortedNames.begin(), SortedNames.end());

    for(int n = 0; n < Ncomp; n++)
    for(int m = 0; m < Ncomp; m++)
    if(Names[n] == SortedNames[m])
    result[m] = n;

    ElementsAreSorted = true;
    for(int n = 0; n < Ncomp; n++)
    if(Names[n] != SortedNames[n])
    ElementsAreSorted = false;

    return result;
}

double ChemicalProperties::GetTracerAmount(int trac)
{
	double TracerAmount = 0.0;
    for(int i = 0; i < Nx; i++)
    for(int j = 0; j < Ny; j++)
    for(int k = 0; k < Nz; k++)
    {
    	TracerAmount += TotalTracer(i,j,k)({trac});
    }
	return TracerAmount;
}

void ChemicalProperties::SetTracerComposition()
{
	//TODO: Error warning for too high tracer composition; calculate total fraction of tracer
    int Bcells = TotalTracer.Bcells();

    for (int x = -Bcells; x < Nx+Bcells; x++)
    for (int y = -Bcells; y < Ny+Bcells; y++)
    for (int z = -Bcells; z < Nz+Bcells; z++)
    for (int trac = 0; trac < NumIsotopes; trac++)
    {
    	if(x != InitialTracerPosition[trac] and x != 500)
    	{
    		TotalTracer(x,y,z)({trac}) = 0.0;
    	}
    	else
    	{
    		TotalTracer(x,y,z)({trac}) = InitialTracerAmount[trac];
    	}
    }
}

void ChemicalProperties::GetInitialTracerAmount()
{
    for(int trac = 0; trac < NumIsotopes ; trac++)
    {
    	minitialTracer[trac] = 0.0;
    }

    for(int i = 0; i < Nx; i++)
    for(int j = 0; j < Ny; j++)
    for(int k = 0; k < Nz; k++)
    for(int trac = 0; trac < NumIsotopes; trac++)
    {
    	minitialTracer[trac] += TotalTracer(i,j,k)({trac});
    }
}

double ChemicalProperties::GetTracerFraction(Composition &Cx, int x, int y, int z, int trac)
{
		double TotalMF = 0.0;
    	for(int alpha = 0; alpha < Nphases; alpha++)
    	{
        	TotalMF += Cx.get_MF(alpha, IsotopeDefined[trac], x,y,z);
    	}

        if(TotalMF > 0.0 and TotalTracer(x,y,z)({trac}) > 0.0)
        return TotalTracer(x,y,z)({trac})/TotalMF;

        else
        return 0.0;
}
}
