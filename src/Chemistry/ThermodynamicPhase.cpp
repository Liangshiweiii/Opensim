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
#include "Chemistry/ThermodynamicPhase.h"
#include "Settings.h"
#include "Compositions.h"

namespace opensim
{
using namespace std;

ThermodynamicPhase::~ThermodynamicPhase(void)
{
    /**Destructor*/
}

ThermodynamicPhase::ThermodynamicPhase(void)
{
    /**Constructor*/
}
    
void ThermodynamicPhase::Initialize(Settings& locSettings,
                         const int boundarysize, PhaseInput& phinput)
{
    Component = phinput.Components;
    Number = phinput.Number;
    vector<string> RefElements = phinput.RefElements;

    /**This function is called from the ChemicalProperties::Initialize() and
    will set the dimension of each storage and fill all parameters, that can
    be taken from the read in input.*/

    thisclassname = "ThermodynamicPhase";
    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    Ncomp = getNcomp();
    Nsubs = getNsubs();
    int ReferenceElementIndex = -1;

    if((int)RefElements.size() > Number)
    {
        for(int comp = 0; comp < Ncomp; comp++)
        if(Component[comp].Name == RefElements[Number])
        {
            Component[comp].Major = true;
            ReferenceElementIndex = Component[comp].Index;
        }

    }

    if(!(ConsIndex.IsAllocated()))
    {
        ConsIndex.Allocate({Nsubs+1});
    }
    else
    {
        ConsIndex.Reallocate({Nsubs+1});
    }
    ConsIndex({0}) = 0;
    for(int s = 0; s < Nsubs; s++)
    {
        Sublattice[s].Initialize();
        ConsIndex({s+1}) = ConsIndex({s}) + Sublattice[s].Ncons;

        /* Check if Reference-Element is present and store position of said
           element, otherwise take first element as reference element of this
           sublattice! */

        int locReferenceElement = 0;
        for(int i = 0; i < Sublattice[s].Ncons; i++)
        if(Sublattice[s].Constituent[i].Index == ReferenceElementIndex)
        {
            locReferenceElement = i;
        }
        Sublattice[s].Reference = locReferenceElement;
    }
    Ncons = getNcons();
    Nsites = getNsites();
    SublIdx = getSublIdx();
    ConsIdx = getConsIdx();
    ConNames = getConNames();
    NsitesWithI = getNsitesWithI();
    dMAdYi = getdMAdYi();
    AnalyticChemicalPotentials = getAnalyticChemicalPotentials();

    if(!(Cmin.IsAllocated()))
    {
        Cmin.Allocate({Ncons});
        Cmax.Allocate({Ncons});
        MolarVolume.Allocate({Ncons});
        Initial.Allocate({Ncons});
        Initial2.Allocate({Ncons});
    }
    else
    {
        Cmin.Reallocate({Ncons});
        Cmax.Reallocate({Ncons});
        MolarVolume.Reallocate({Ncons});
        Initial.Reallocate({Ncons});
        Initial2.Reallocate({Ncons});
    }
    //Composition.Allocate(Nx, Ny, Nz, {Ncons}, boundarysize);
    //MoleFractions.Allocate(Nx, Ny, Nz, {Ncomp}, boundarysize);


    //SublDot.resize(Nsubs);
    /*for(int sub = 0; sub < Nsubs; sub++)
    {
        int subNcons = Sublattice[sub].Ncons;
        //SublDot[sub].Allocate(Nx,Ny,Nz,{subNcons,subNcons,3,3,3},boundarysize);
    }*/
    SiteFractions.resize(Ncons);
    for(int con = 0; con < Ncons; con++)
    {
        SiteFractions[con].Min = 0.0;
        SiteFractions[con].Max = 1.0;
    }

    if(locSettings.CP.diffMod == ChemicalProperties::DiffusionModel::EquilibriumPartitioning)
    {
        for(int comp = 0; comp < Ncomp; comp++)
        {
            Component[comp].isStoichiometric = false;
        }
    }
    else if(locSettings.CP.diffMod == ChemicalProperties::DiffusionModel::FiniteInterfaceDissipation)
    {
        for(int comp = 0; comp < Ncomp; comp++)
        {
            Component[comp].isStoichiometric
            = isStoichiometric(Component[comp].Index);
        }
    }

    //InitialTracerPhaseComposition.Allocate({Ncons});//TODO: DELETE!
}

Tensor<double, 1> ThermodynamicPhase::Site2MoleFrac(Composition& Cx, int i, int j, int k)
{
    //TODO: Please use Composition::Site2MoleFrac
    /**This tensor */
    Tensor<double, 1> tempFrac;
    tempFrac.Allocate({Ncomp});
    for(int n = 0; n < Ncomp; n++)
    {
        double mfrac = 0.0;
        double div = 0.0;
        int cons = 0;
        for(int s = 0; s < Nsubs; s++)
        {
            double Site = Sublattice[s].Site;
            double vac = 0.0;
            for(int nc = 0; nc < Sublattice[s].Ncons; nc++)
            {

                if(Sublattice[s].Constituent[nc].Index == Component[n].Index
                   and !(Sublattice[s].Constituent[nc].isVacancy))
                {
                    mfrac += Cx.get_CF(Number,cons,i,j,k)*Site;
                }
                else if(Sublattice[s].Constituent[nc].isVacancy)
                {
                    vac += Cx.get_CF(Number,cons,i,j,k);
                }
                cons++;
            }
            div += (1.0-vac)*Site;
        }
        tempFrac[n] = mfrac/div;
    }
    return tempFrac;
}

void ThermodynamicPhase::CalculateMoleFractions(Composition& Cx)
{
    /**This function will loop over the whole simulation domain and calculate
    and update the mole fraction in each point from the stored site fractions.*/
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Cx.Phase,0,)
    {
        Cx.set_MFs(Number,Site2MoleFrac(Cx, i, j, k),i,j,k);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

int ThermodynamicPhase::getNcons(void)
{
    /**This function will return the total number of constituents for this
    phase.*/
    int tmp = 0;
    for(int s = 0; s < Nsubs; s++)
    {
        tmp += Sublattice[s].Ncons;
    }
    return tmp;    
}

int ThermodynamicPhase::getNsubs(void)
{
    /**This function will return the number of sublattices for this phase.*/
    return Sublattice.size();    
}

int ThermodynamicPhase::getNcomp(void)
{
    /** This function returns the number of components (!not in the present
    phase, but in the total system) as there is no access to
    ChemicalProperties.Ncomp from within this class.*/
    return Component.size();    
}

double ThermodynamicPhase::getNsites(void)
{
    /** This function returns the sum of site coefficients of all
    sublattices. For the phase (A,B)_1 (A,B,C)_3 (C)_5 this function will
    return 1 + 3 + 5 = 9.*/
    double locNsites = 0.0;
    for(int n = 0; n < Nsubs; n++)
    {
        locNsites += Sublattice[n].Site;
    }
    return locNsites;
}

double ThermodynamicPhase::Nmoles(Composition& Cx, int x, int y, int z)
{
    /** This function returns the number of moles in the present phase (not
    weighted with the phase field). Does not return the overall number of
    moles over all phases in the point [x,y,z].*/
    int counter = 0;
    double tempNmoles = 0.0;
    for(int n = 0; n < Nsubs; n++)
    {
        double tempSite = Sublattice[n].Site;
        double temp = 0.0;
        for(int i = 0; i < Sublattice[n].Ncons; i++)
        {
            if(Sublattice[n].Constituent[i].Index >= 0)
            {
                temp += Cx.get_CF(Number,counter,x,y,z);
            }
            counter++;
        }
        tempNmoles += temp*tempSite;
    }
    return tempNmoles;
}

std::vector<std::string> ThermodynamicPhase::getConNames(void)
{
    /** This vector of strings holds all the Names of each constituent in
    this phase. Can be used for printouts. */
    std::vector<std::string> tempNames;
    for(int sub = 0; sub < Nsubs; sub++)
    for(int con = 0; con < Sublattice[sub].Ncons; con++)
    {
        tempNames.push_back(Sublattice[sub].Constituent[con].Name);
    }
    return tempNames;
}

std::vector<int> ThermodynamicPhase::getSublIdx(void)
{
    /** For the example sublattice (A,B,C)(D,E)(F)(G,H), SublIdx will hold
    the values [0,3,5,6,8]. This helps when looping over all constituents of
    a phase, but where the sublattice has to be identified. Constituent 0-2,
    3-4, 5, 6-7 all share the same sublattice. Use as:
            for(int sub = 0; sub < Nsubs; sub++)
            for(int con = SublIdx[sub]; con < SublIdx[sub+1]; con++)*/
    std::vector<int> tempSublIdx(Nsubs+1);
    tempSublIdx[0] = 0;
    for(int sub = 0; sub < Nsubs; sub++)
    {
        tempSublIdx[sub+1] = tempSublIdx[sub]+Sublattice[sub].Ncons;
    }
    return tempSublIdx;
}

std::vector<int> ThermodynamicPhase::getConsIdx(void)
{
    /**This function returns a vector of the size of all Constituents of the
    phase, that holds the Index of each Constituent to identify each
    Constituent as a Component. If you want to know the Index of one
    Constituent 'con', you have to use ConsIdx[con].*/
    std::vector<int> tempConsIdx(Ncons,-1);
    for(int com = 0; com < Ncomp; com++)
    {
        int counter = 0;
        for(int sub = 0; sub < Nsubs; sub++)
        for(int con = 0; con < Sublattice[sub].Ncons; con++)
        {
            if(Component[com].Index == Sublattice[sub].Constituent[con].Index)
            {
                tempConsIdx[counter] = Component[com].Index;
            }
            counter++;
        }
    }
    return tempConsIdx;
}

int ThermodynamicPhase::Idx2Cons(int sub, int Idx)
{
    /**This function returns the constituent number of the given component index
    Idx*/
    int result = -1;
    for(int con = SublIdx[sub]; con < SublIdx[sub+1]; con++)
    if(ConsIdx[con] == Idx)
    result = con;
    if(result < 0)
    {
        cout << "Error in Idx2Cons("<<sub<<","<<Idx<<")" << endl;
        //exit(1);
    }
    return result;
}

std::vector<double> ThermodynamicPhase::yTotal(Composition& Cx, int x, int y, int z)
{
    /** This function sums up over the product site*sitefraction for each
    sublattice for each component at the point [x,y,z]. If you divide yTotal
    my the number of moles, you will get the mole fractions at this point.*/
    std::vector<double> temp(Ncomp+1, 0.0);
    std::vector<int> locConsIdx = ConsIdx;
    for(int sub = 0; sub < Nsubs; sub++)
    for(int con = SublIdx[sub]; con < SublIdx[sub+1]; con++)
    {
        if(locConsIdx[con] >= 0)
        {
            temp[locConsIdx[con]]  += Sublattice[sub].Site
                                      *Cx.get_CF(Number,con,x,y,z);
        }
        else
        {
            temp[Ncomp]  += Sublattice[sub].Site*Cx.get_CF(Number,con,x,y,z);
        }
    }
    return temp;
}

std::vector<double> ThermodynamicPhase::getNsitesWithI(void)
{
    /**For the example phase (A,B)_1 (A)_3 (A,B)_5 this function will
    hold the following values: NsitesWithI[A] = 9, NsitesWithI[B] = 6,
    NsitesWithI[C] = 0. For the total summation of all Sites, see Nsites.*/
    std::vector<double> temp(Ncomp+1, 0.0);
    std::vector<int> locConsIdx = ConsIdx;
    for(int sub = 0; sub < Nsubs; sub++)
    for(int con = SublIdx[sub]; con < SublIdx[sub+1]; con++)
    {
        if(locConsIdx[con] >= 0)
        {
            temp[locConsIdx[con]]  += Sublattice[sub].Site;
        }
        else
        {
            temp[Ncomp]  += Sublattice[sub].Site;
        }
    }
    return temp;
}


std::vector<double> ThermodynamicPhase::unityWithoutVA(Composition& Cx, int x, int y, int z)
{
    /**This function returns the summation of all site fractions of elements
    in one sublattice, not including vacant sites. This will give the
    correct y_Va for each sublattice in one vector, .*/
    int locNsubs = Nsubs;
    std::vector<int> locConsIdx = ConsIdx;
    std::vector<double> temp(locNsubs, 0.0);
    for(int sub = 0; sub < locNsubs; sub++)
    for(int con = SublIdx[sub]; con < SublIdx[sub+1]; con++)
    if(locConsIdx[con] >= 0)
    {
        temp[sub] += Cx.get_CF(Number,con,x,y,z);
    }
    return temp;
}

double ThermodynamicPhase::Nm(Composition& Cx, int x, int y, int z)
{
    /**might be obsolete*/
    double tempSites = 0.0;
    for(int n = 0; n < Nsubs; n++)
    {
        tempSites += Sublattice[n].Site;
    }
    double tempNm = tempSites / Nmoles(Cx,x,y,z);
    return tempNm;
}

bool ThermodynamicPhase::isInterstitial(int comp)
{
    /**This boolean value will identify if the component is an interstitial
    component (true, if present on sublattice with vacancies) or an
    substitutional component (false, if not present on any sublattice with
    vacancies).*/
    bool temp = false;
    for(int sub = 0; sub < Nsubs; sub++)
    if(Sublattice[sub].hasVacancies
       and Sublattice[sub].isElementPresent(Component[comp].Index))
    {
        temp = true;
    }
    return temp;
}

bool ThermodynamicPhase::isEndmember(int CIdxI)
{
    /**This boolean value will identify if the component I is an endmember in
    this phase. That means if component I exists on every sublattice, which
    doesn't have vacancies. CIdxI has to be the index of the component!!*/

    bool Iispresent = false;
    bool Iisendmember = false;

    for(int sub = 0; sub < Nsubs; sub++)
    if(Sublattice[sub].isElementPresent(CIdxI))
    Iispresent = true;

    if(Iispresent)
    {
        Iisendmember = true;
        for(int sub = 0; sub < Nsubs; sub++)
        if((not Sublattice[sub].isElementPresent(CIdxI))
        and(not Sublattice[sub].hasVacancies)
        /*and(Sublattice[sub].Ncons > 1)*/)
        {// no I nor Va on sublattice -> I not an endmember (if Ncons > 1)
            Iisendmember = false;
        }
    }
    return Iisendmember;
}

bool ThermodynamicPhase::ispairwithVA(int CIdxI)
{
    /**This boolean value will identify if the component I shares its sublattice
    with vacancy VA on at least one sublattice. CIdxI has to be the index of the
    component!!*/

    bool IshareswithVA = false;

    for(int sub = 0; sub < Nsubs; sub++)
    if((Sublattice[sub].isElementPresent(CIdxI)) and Sublattice[sub].isElementPresent(-1))
    IshareswithVA = true;

    return IshareswithVA;
}

bool ThermodynamicPhase::isStoichiometric(int CIdxI)
{
    /**This boolean value will identify if the component I is a stoichiometric
    component, as it shares its sublattice without any other component.*/

    bool aloneinsubl = true;
    bool nointerstitials = true;

    bool isStoich = false;

    for(int sub = 0; sub < Nsubs; sub++)
    {
        if((Sublattice[sub].isElementPresent(CIdxI)) and Sublattice[sub].Ncons > 1)
        {
            aloneinsubl = false;
        }

        if(Sublattice[sub].hasVacancies and Sublattice[sub].Ncons > 1)
        {
            nointerstitials = false;
        }
    }

    if(aloneinsubl and nointerstitials)
    isStoich = true;

    return isStoich;
}

bool ThermodynamicPhase::get_isStoichiometric(void)
{
    /***/

    bool isStoich = false;

    for(int n = 0; n < Ncons; n++)
    if(ConsIdx[n] > -1 and isStoichiometric(ConsIdx[n]))
    {
        isStoich = true;
    }

    return isStoich;
}

vector<vector<double> > ThermodynamicPhase::getdMAdYi(void)
{
    /**This matrix is populated with values for dM_A/dY_i, which are necessary
    for calculation of chemical potentials*/

    vector<vector<double> > result;

    for(int comp = 0; comp < Ncomp; comp++)
    {
        vector<double> tempcons(Ncons, 0.0);

        int counter = 0;
        for(int sub = 0; sub < Nsubs; sub++)
        for(int con = 0; con < Sublattice[sub].Ncons; con++)
        {
            if(Sublattice[sub].Constituent[con].Index == Component[comp].Index)
            {
                tempcons[counter] = Sublattice[sub].Site;
            }
            counter++;
        }

        result.push_back(tempcons);
    }

    return result;
}

bool ThermodynamicPhase::getAnalyticChemicalPotentials(void)
{
    /**This boolean determines, whether chemical potentials can be calculated
    analytically or only by numeric calculation*/

    bool result = true;

    for(int A = 0; A < Ncomp; A++)
    {
        int CIdxI = Component[A].Index;
        if(!(isEndmember(CIdxI) or ispairwithVA(CIdxI)))
        result = false;

    }

    //check if all components are endmembers or share a sublattice with VA

    return result;
}

}
