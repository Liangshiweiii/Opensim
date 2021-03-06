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

#include "Nucleation.h"
#include "Settings.h"
#include "PhaseField.h"
#include "Temperatures.h"
#include "DrivingForce.h"
#include "InterfaceEnergy.h"
#include "Tools/UserInterface.h"
#include "Orientations.h"
#include "Info.h"
#include "Mechanics/Storages/SymmetryVariants.h"

namespace opensim
{
using namespace std;

Nucleation::Nucleation(const Settings& locSettings, const std::string InputFileName)
{
    this->Initialize(locSettings);

    if(InputFileName == "NONE")
    {
        this->ReadInput();
    }
    else
    {
        this->ReadInput(InputFileName);
    }
}

void Nucleation::Initialize(const Settings& locSettings)
{
    thisclassname = "Nucleation";
    //DefaultInputFileName = ProjectInputDir + "NucleationInput.opi";

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;
    Nphases = locSettings.Nphases;

    NucleationParameters.Allocate(Nphases, Nphases);
    GeneratedParticles.Allocate(Nphases, Nphases);

    SeedX = 4556;
    SeedY = 21;
    SeedZ = 359;

    ofstream NucFile;
    NucFile.open ("NucleationStatistics.dat", ios::out);
    NucFile << setw(14) << "TimeStep"
            << setw(8)  << "PFindex"
            << setw(8)  << "Nphase"
            << setw(8)  << "Mphase"
            << setw(6)  << "x"
            << setw(6)  << "y"
            << setw(6)  << "z"
            << setw(12) << "Q1"
            << setw(12) << "Q2"
            << setw(12) << "Q3"
            << setw(12) << "dGnuc"
            << setw(12) << "dGmin"
            << setw(10) << "Status" << endl;

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void Nucleation::ReadInput(const std::string InputFileName)
{
    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        Info::WriteExit("File " + InputFileName + " could not be opened",
                        thisclassname, "ReadInput()");
        exit(16);
    }

    Info::WriteBlankLine();
    Info::WriteLineInsert("Nucleation");
    Info::WriteStandard("Source", InputFileName);
    
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
    // Reading Nucleation Allowance Booleans and corresponding nucleation densities:
    for(int n = 0; n < Nphases; n++)
    for(int m = 0; m < Nphases; m++)
    {
        stringstream converterN;
        converterN << n;
        stringstream converterM;
        converterM << m;
        string counterN = converterN.str();
        string counterM = converterM.str();
        
        string counter = counterN + "_" + counterM;
        
        string userinput = UserInterface::ReadParameterS(inp, moduleLocation, string("Allowed_")
                                                              + counter);

        std::transform(userinput.begin(), userinput.end(),
                       userinput.begin(), ::toupper);

        //Read user defined nucleation modes for phase M in phase N
        // "YES"  allows nucleation globally in the simulation box
        // "NO"   forbids nucleation globally in the simulation box
        // "GB"   only allows nucleation on grain boundaries where M is present
        // "BULK" only allows nucleation in bulk regions of phase

        if(userinput != "NO")
        {
            NucleationParameters(n, m).Allowed = 1;
            NucleationParameters(n, m).Density        = UserInterface::ReadParameterD(inp, moduleLocation, string("NucDens_") + counter);

            // Reading Particle Distribution
            NucleationParameters(n, m).DistMu         = UserInterface::ReadParameterD(inp, moduleLocation, string("DistMu_") + counter);
            NucleationParameters(n, m).DistSigma      = UserInterface::ReadParameterD(inp, moduleLocation, string("DistSig_") + counter);
            NucleationParameters(n, m).Tmin           = UserInterface::ReadParameterD(inp, moduleLocation, string("Tmin_") + counter);
            NucleationParameters(n, m).Tmax           = UserInterface::ReadParameterD(inp, moduleLocation, string("Tmax_") + counter);
            NucleationParameters(n, m).Generated      = false;

            // Reading Particle Orientation
            string tmpString = UserInterface::ReadParameterS(inp, moduleLocation, string("Orient_") + counter);
            if(tmpString == "Random")
            {
                NucleationParameters(n, m).OrientationMode = Random;
            }
            else
            {
                NucleationParameters(n, m).OrientationMode = 0;
            }
            
            tmpString = UserInterface::ReadParameterS(inp, moduleLocation, string("Variant_") + counter);            
            if(tmpString == "Single")
            {
                NucleationParameters(n, m).Variants = 0;
            }
            else
            {
                NucleationParameters(n, m).Variants = 1;
            }
            // Reading Nucleation Sites
            if(userinput == "YES")
            {
                NucleationParameters(n, m).Mode = 0;
            }
            else if(userinput == "BULK")
            {
                NucleationParameters(n, m).Mode = 1;
            }
            else if(userinput == "GB")
            {
                NucleationParameters(n, m).Mode = 2;
            }
        }
        else if(userinput == "NO")
        {
            NucleationParameters(n, m).OrientationMode = 0;
            NucleationParameters(n, m).Density         = 0.0;
            NucleationParameters(n, m).Generated       = false;
        }
        else
        {
            Info::WriteExit("Nucleation mode for phase pair " + counterN + "_"
                            + counterM + " could not be read", thisclassname,
                            "ReadInput()");
            exit(16);
        }
    }
    
    inp.close();

    Info::WriteLine();
}

void Nucleation::GenerateRandomSeeds(void)
{
    int min = 1;
    int max = 10000;

    SeedX = min + (rand() % (int)(max - min + 1));
    SeedY = min + (rand() % (int)(max - min + 1));
    SeedZ = min + (rand() % (int)(max - min + 1));
}

void Nucleation::Clear()
{
    for(int n = 0; n < Nphases; n++)
    for(int m = 0; m < Nphases; m++)
    {
        GeneratedParticles(n,m).clear();
    }
}

vector<vector<int> > Nucleation::GetNucleationSites(PhaseField& Phi)
{
    /** This function return a vector of a coordinates-vector for all points
    where nucleation events have happend. Condition is planted = true and
    Fractions < 1E-1*/

    vector<vector<int> > result;

    for(int n = 0; n < Nphases; n++)
    for(int m = 0; m < Nphases; m++)
    for(auto ind = GeneratedParticles(n, m).begin();
             ind != GeneratedParticles(n, m).end(); ind++)
    {
        vector<int> temp(3, 0);
        temp[0] = ind->x;
        temp[1] = ind->y;
        temp[2] = ind->z;

        if(ind->planted)
        {
            int Phase = Phi.FieldsStatistics[ind->PFindex].Phase;
            double Fractions = Phi.Fractions(temp[0],temp[1],temp[2])({Phase});

            if(Fractions < 1E-1)
            {
                result.push_back(temp);
            }
        }
    }

    return result;
}


bool Nucleation::IsShielded(std::vector<NucSite> LocGrainStorage, const int i, const int j, const int k, const double shielding) const
{
    for (auto it = LocGrainStorage.begin(); it != LocGrainStorage.end(); ++it)
    {
        const double xdis = std::min(std::fabs(i - it->x), std::min( std::fabs(i - it->x + Nx), std::fabs(i - it->x - Nx)));
        const double ydis = std::min(std::fabs(j - it->y), std::min( std::fabs(j - it->y + Ny), std::fabs(j - it->y - Ny)));
        const double zdis = std::min(std::fabs(k - it->z), std::min( std::fabs(k - it->z + Nz), std::fabs(k - it->z - Nz)));
        const double distanceSquare = xdis*xdis + ydis*ydis + zdis*zdis;
        if(distanceSquare < shielding*shielding) return true;
    }
    return false;
}

void Nucleation::GenerateNucleationSites(PhaseField& Phase, Temperature& Tx)
{
    // Using mersenne-twister 64 bit pseudo-random generator engine:
    mt19937_64 SizeGenerator(53);

    mt19937_64 xPosGenerator(SeedX);
    mt19937_64 yPosGenerator(SeedY);
    mt19937_64 zPosGenerator(SeedZ);

    mt19937_64 OrientGenerator1(45);
    mt19937_64 OrientGenerator2(697);
    mt19937_64 OrientGenerator3(255);

    for(int n = 0; n < Nphases; n++)
    for(int m = 0; m < Nphases; m++)
    if(NucleationParameters(n, m).Allowed and
       NucleationParameters(n, m).Tmin < Tx(0, 0, 0) and
       NucleationParameters(n, m).Tmax > Tx(0, 0, 0) and
       !NucleationParameters(n, m).Generated)
    {
        double DistMu = NucleationParameters(n, m).DistMu;
        double DistSigma = NucleationParameters(n, m).DistSigma;

        normal_distribution <double> SizeDistribution(DistMu, DistSigma);

        uniform_int_distribution <int> xPosDistribution(0, Nx - 1);
        uniform_int_distribution <int> yPosDistribution(0, Ny - 1);
        uniform_int_distribution <int> zPosDistribution(0, Nz - 1);

        uniform_int_distribution <int> Q1Distribution(0, 180);
        uniform_int_distribution <int> Q2Distribution(0, 180);
        uniform_int_distribution <int> Q3Distribution(0, 180);

        vector<double> PhaseFractions(Phase.Nphases,0.0);
        double TotalVolume = Phase.Nx*Phase.Ny*Phase.Nz;
        for(unsigned int idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
        {
            PhaseFractions[Phase.FieldsStatistics[idx].Phase] += Phase.FieldsStatistics[idx].Volume/TotalVolume;
        }

        double dNsites = NucleationParameters(n, m).Density*
                         PhaseFractions[m]*
                         Nx*Ny*Nz*
                         Phase.dx*Phase.dx*Phase.dx;
                        
        int Nsites = dNsites;
        
        if(Nsites)
        {
            int    part = 0;
            int    attempts = 0;
            double shielding = 2.0*Phase.iWidth;
            while (part < Nsites and attempts < Nx*Ny*Nz)
            {
                double weight = SizeDistribution(SizeGenerator);

                int Xpos = xPosDistribution(xPosGenerator);
                int Ypos = yPosDistribution(yPosGenerator);
                int Zpos = zPosDistribution(zPosGenerator);

                bool locmode = false;

                switch (NucleationParameters(n, m).Mode)
                {
                    case 0:
                    {
                        //Global nucleation
                    	locmode = true;
                        break;
                    }
                    case 1:
                    {
                        //Bulk nucleation
                        if(!(Phase.Fields(Xpos,Ypos,Zpos).flag))
                        	locmode = true;
                        break;
                    }
                    case 2:
                    {
                        //Grain boundary nucleation
                        if(Phase.Fields(Xpos,Ypos,Zpos).flag)
                        	locmode = true;
                        break;
                    }                    
                    default:
                    {
                        //Fail safe
                        break;
                    }
                }

                if (!IsShielded(GeneratedParticles(n, m), Xpos, Ypos, Zpos, shielding) and (weight > 0) and 
                    (Phase.Fractions(Xpos, Ypos, Zpos)({m}) > 0.5) and locmode)
                {
                    NucSite locSite;
                    locSite.x = Xpos;
                    locSite.y = Ypos;
                    locSite.z = Zpos;
                    locSite.radius = weight;

                    if(NucleationParameters(n, m).OrientationMode == Random)
                    {
                        locSite.Q1 = Q1Distribution(OrientGenerator1) * Pi/180.0;
                        locSite.Q2 = Q2Distribution(OrientGenerator2) * Pi/180.0;
                        locSite.Q3 = Q3Distribution(OrientGenerator3) * Pi/180.0;
                    }
                    else
                    {
                        locSite.Q1 = 0.0;
                        locSite.Q2 = 0.0;
                        locSite.Q3 = 0.0;
                    }
                    std::stringstream message;
                    message << "Nucleation: Generated seed particle " << part << " at [" 
                            << Xpos << ", " << Ypos << ", " << Zpos << "] and effective radius of " << weight;
                    Info::WriteSimple(message.str());
                    GeneratedParticles(n, m).push_back(locSite);
                    part ++;
                }
                else
                {
                   attempts++;
                }
            }
            NucleationParameters(n, m).Generated = true;
            std::stringstream message2;
            message2 << "Nucleation: Generated " << part << " nucleation sites for phase " << n << " in phase " << m << ".";
            Info::WriteSimple(message2.str());
        }
        else
        {
            std::stringstream message3;
            message3 << "Nucleation: Too low particles density! No nucleation sites were generated for phase " 
                     << n << " in phase " << m << ".";
            Info::WriteSimple(message3.str());
        }
    }
}

void Nucleation::WriteStatistics(int tStep, int PFindex, int NuleationPhase, 
                                    int MatrixPhase, int Xpos, int Ypos, int Zpos, 
                                    double Q1, double Q2, double Q3, double dGnuc,
                                    double dGmin, string status) const
{
    ofstream NucFile;
    NucFile.open ("NucleationStatistics.dat", ios::app);
    
    NucFile << setw(14)  << setprecision(6) << tStep
            << setw( 8)  << setprecision(6) << PFindex
            << setw( 8)  << setprecision(6) << NuleationPhase            
            << setw( 8)  << setprecision(6) << MatrixPhase
            << setw( 6)  << setprecision(5) << Xpos
            << setw( 6)  << setprecision(5) << Ypos
            << setw( 6)  << setprecision(5) << Zpos
            << setw(12)  << setprecision(4) << Q1
            << setw(12)  << setprecision(4) << Q2
            << setw(12)  << setprecision(4) << Q3            
            << setw(12)  << setprecision(6) << dGnuc
            << setw(12)  << setprecision(6) << dGmin
            << setw(10)  << status          << endl;
    NucFile.close();
}

void Nucleation::PlantNuclei(PhaseField& Phase, int tStep)
{    
    for(int n = 0; n < Nphases; n++)
    for(int m = 0; m < Nphases; m++)
    for(auto ind = GeneratedParticles(n, m).begin();
             ind != GeneratedParticles(n, m).end(); ind++)
    {
        if(Phase.Fractions(ind->x, ind->y, ind->z)({n}) == 0.0 and 
           Phase.Fractions(ind->x, ind->y, ind->z)({m}) >= 0.3)
        {
            int locIndex = Phase.PlantGrainNucleus(n, ind->x, ind->y, ind->z);
            EulerAngles locAngles({ind->Q1, ind->Q2, ind->Q3}, XYZ);
            Phase.FieldsStatistics[locIndex].Orientation = locAngles.getQuaternion();
            ind->PFindex = locIndex;
            ind->planted = true;
        }
        else 
        {
            WriteStatistics(tStep, -1, n, m, ind->x, ind->y, ind->z, ind->Q1, ind->Q2, ind->Q3, 0, 0, "Deleted");
            GeneratedParticles(n, m).erase(ind);
            ind--;
        }
    }
}

void Nucleation::PlantNuclei(PhaseField& Phase, SymmetryVariants& SV, int tStep)
{
    // Using mersenne-twister 64 bit pseudo-random generator engine:
    //mt19937_64 VariantsGenerator(tStep);
    
    for(int n = 0; n < Nphases; n++)
    for(int m = 0; m < Nphases; m++)
    for(auto ind = GeneratedParticles(n, m).begin();
             ind != GeneratedParticles(n, m).end(); )
    {
        //uniform_int_distribution <int> VariantSelector(0, SV.Nvariants(n)-1);
        
        if(Phase.Fractions(ind->x, ind->y, ind->z)({n}) == 0.0 and 
           Phase.Fractions(ind->x, ind->y, ind->z)({m}) >= 0.3)
        {
            //int v = VariantSelector(VariantsGenerator);
            for(int v = 0; v < SV.Nvariants(n); v++)
            //if(VariantProbability(VariantsSelector) > 0.5)
            {
                int locIndex = Phase.PlantGrainNucleus(n, ind->x, ind->y, ind->z);
                EulerAngles locAngles({ind->Q1, ind->Q2, ind->Q3}, XYZ);
                Phase.FieldsStatistics[locIndex].Orientation = locAngles.getQuaternion();
                Phase.FieldsStatistics[locIndex].Variant = v;
                ind->PFindices.push_back(locIndex);
                //ind->PFindex = locIndex;
            }
            ind->planted = true;
            ind++;
        }
        else 
        {
            WriteStatistics(tStep, -1, n, m, ind->x, ind->y, ind->z, ind->Q1, ind->Q2, ind->Q3, 0, 0, "Deleted");
            ind = GeneratedParticles(n, m).erase(ind);
        }
    }
}

void Nucleation::PlantNucleiGB(PhaseField& Phase, SymmetryVariants& SV, int tStep)
{
    // Using mersenne-twister 64 bit pseudo-random generator engine:
    //mt19937_64 VariantsGenerator(tStep);
    
    for(int n = 0; n < Nphases; n++)
    for(int m = 0; m < Nphases; m++)
    for(auto ind = GeneratedParticles(n, m).begin();
             ind != GeneratedParticles(n, m).end(); )
    {
        //uniform_int_distribution <int> VariantSelector(0, SV.Nvariants(n)-1);
        
        if(Phase.Fractions(ind->x, ind->y, ind->z)({n}) >= 0.3 and 
           Phase.Fractions(ind->x, ind->y, ind->z)({m}) >= 0.3)
        {
            //int v = VariantSelector(VariantsGenerator);
            for(int v = 0; v < SV.Nvariants(n); v++)
            //if(VariantProbability(VariantsSelector) > 0.5)
            {
                int locIndex = Phase.PlantGrainNucleus(n, ind->x, ind->y, ind->z);
                EulerAngles locAngles({ind->Q1, ind->Q2, ind->Q3}, XYZ);
                Phase.FieldsStatistics[locIndex].Orientation = locAngles.getQuaternion();
                Phase.FieldsStatistics[locIndex].Variant = v;
                ind->PFindices.push_back(locIndex);
                //ind->PFindex = locIndex;
            }
            ind->planted = true;
            ind++;
        }
        else if (Phase.Fractions(ind->x, ind->y, ind->z)({n}) > 0.7)
        {
            WriteStatistics(tStep, -1, n, m, ind->x, ind->y, ind->z, ind->Q1, ind->Q2, ind->Q3, 0, 0, "Deleted");
            ind = GeneratedParticles(n, m).erase(ind);
        }
        else
        {
            ind++;
            ind->planted = false;
        }
    }
}

void Nucleation::CheckNuclei(PhaseField& Phase, InterfaceEnergy& IE, DrivingForce& dG, int tStep)
{
    for(int n = 0; n < Nphases; n++)
    for(int m = 0; m < Nphases; m++)
    for(auto ind = GeneratedParticles(n, m).begin();
             ind != GeneratedParticles(n, m).end(); )
    {
    if(ind->planted)
    {
        double dGmin = 2.0*IE.MinIntEnergy(n, m)/ind->radius;
        double dGloc = 0.0;

        for(auto it = dG.Raw(ind->x, ind->y, ind->z).begin(); 
                 it != dG.Raw(ind->x, ind->y, ind->z).end(); ++it)
        {
            if(it->indexA == ind->PFindex)
            {
                dGloc += it->value;
            }
            if(it->indexB == ind->PFindex)
            {
                dGloc -= it->value;
            }
        }
        //cout << dGloc << " " << dGmin << endl;
        if(dGloc < dGmin)
        {
            ind->planted = false;
            Phase.FieldsStatistics[ind->PFindex].Exist = false;
            for(auto it = dG.Raw(ind->x, ind->y, ind->z).begin();
                     it != dG.Raw(ind->x, ind->y, ind->z).end(); )
            {
                if(it->indexA == ind->PFindex or it->indexB == ind->PFindex)
                {
                    it = dG.Raw(ind->x, ind->y, ind->z).erase(it);
                }
                else
                {
                    ++it;
                }
            }
            ++ind;
        }
        else
        {
            WriteStatistics(tStep, ind->PFindex, n, m, ind->x, ind->y, ind->z, ind->Q1, ind->Q2, ind->Q3, dGloc, dGmin, "Planted");                
            ind = GeneratedParticles(n, m).erase(ind);
        }
    }
    else
    {
        //TODO: check for consistency (this prevents freeze of simulation)
        break;
    }
    }
}

void Nucleation::CheckNuclei(PhaseField& Phase, InterfaceEnergy& IE, DrivingForce& dG, SymmetryVariants& SV, int tStep)
{
    // Using mersenne-twister 64 bit pseudo-random generator engine:
    mt19937_64 VariantsGenerator(tStep);
    
    for(int n = 0; n < Nphases; n++)
    for(int m = 0; m < Nphases; m++)
    for(auto ind = GeneratedParticles(n, m).begin();
             ind != GeneratedParticles(n, m).end(); )
    if(ind->planted)
    {
        double dGmin = 2.0*IE.MinIntEnergy(n, m)/ind->radius;
        vector<double> dGloc(SV.Nvariants(n), 0.0);
        
        for(auto it = dG.Raw(ind->x, ind->y, ind->z).begin(); 
                 it != dG.Raw(ind->x, ind->y, ind->z).end(); ++it)
        {
            for(int v = 0; v < SV.Nvariants(n); v++)
            {
                if(it->indexA == ind->PFindices[v])
                {
                    dGloc[v] += it->value;
                }
                if(it->indexB == ind->PFindices[v])
                {
                    dGloc[v] -= it->value;
                }
            }
        }
        
        vector<double> Probability(SV.Nvariants(n), 0);
        double minProbability = 0.0;
        double maxProbability = 1.0;
        uniform_real_distribution <double> VariantSelector(minProbability, maxProbability);
        for(int v = 0; v < SV.Nvariants(n); v++)
        {
            if(dGloc[v] > dGmin)
            {
                Probability[v] = VariantSelector(VariantsGenerator);
                //WriteStatistics(tStep, ind->PFindices[v], n, m, ind->x, ind->y, ind->z, ind->Q1, ind->Q2, ind->Q3, dGloc[v], dGmin, "Planted");                
            }
        }
        int mostProbableVariant = 0;
        for(int v = 0; v < SV.Nvariants(n); v++)
        {
            if(Probability[v] > minProbability)
            {
                minProbability = Probability[v];
                mostProbableVariant = v;
            }
        }
        
        int NucleiCounter = SV.Nvariants(n);
        for(int v = 0; v < SV.Nvariants(n); v++)
        {
            if(v != mostProbableVariant)
            {
                NucleiCounter--;
                
                Phase.FieldsStatistics[ind->PFindices[v]].Exist = false;
                Phase.FieldsStatistics[ind->PFindices[v]].Stage = 0;
                
                for(auto it = dG.Raw(ind->x, ind->y, ind->z).begin();
                         it != dG.Raw(ind->x, ind->y, ind->z).end(); )
                {
                    if(it->indexA == ind->PFindices[v] or 
                       it->indexB == ind->PFindices[v])
                    {
                        it = dG.Raw(ind->x, ind->y, ind->z).erase(it);
                    }
                    else
                    {
                        ++it;
                    }
                }
            }
            else if (Probability[v] != 0.0)
            {
                WriteStatistics(tStep, ind->PFindices[v], n, m, ind->x, ind->y, ind->z, ind->Q1, ind->Q2, ind->Q3, dGloc[v], dGmin, "Planted");                
            }
        }
        
        if(NucleiCounter != 0)
        {
            ind = GeneratedParticles(n, m).erase(ind);
        }
        else
        {            
            ind->planted = false;
            ind->PFindices.resize(0);
            ind++;
        }
    }
    else
    {
        ind++;
    }
}

void Nucleation::Write(int tStep)
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"Nucleation_", tStep, ".dat");
    fstream out(FileName.c_str(), ios::out);

    if (!out)
    {
        Info::WriteExit("File/" + FileName + "could not be created", thisclassname);
        exit(1);
    }

    for(int n = 0; n < Nphases; ++n)
    {
        for(int m = 0; m < Nphases; ++m)
        {
            out << GeneratedParticles(n,m).size() << " " ;
        }
        out << endl;
    }
    for(int n = 0; n < Nphases; ++n)
    for(int m = 0; m < Nphases; ++m)
    for(unsigned int i = 0; i < GeneratedParticles(n,m).size(); ++i)
    {
        out << GeneratedParticles(n,m)[i].x  << " " 
            << GeneratedParticles(n,m)[i].y  << " " 
            << GeneratedParticles(n,m)[i].z  << " " << endl;
        out << GeneratedParticles(n,m)[i].Q1 << " " 
            << GeneratedParticles(n,m)[i].Q2 << " " 
            << GeneratedParticles(n,m)[i].Q3 << " " << endl;
        out << GeneratedParticles(n,m)[i].radius << endl;
    }
    out.close();
}

void Nucleation::Read(int tStep)
{
    string FileName = UserInterface::MakeFileName(RawDataDir,"Nucleation_", tStep, ".dat");
    fstream inp(FileName.c_str(), ios::in);

    if (!inp)
    {
        Info::WriteExit("File/" + FileName + "could not be opened", thisclassname);
        exit(1);
    };

    for(int n = 0; n < Nphases; ++n)
    for(int m = 0; m < Nphases; ++m)
    {
        int tmp;
        inp >> tmp;

        if(tmp)
        {
            GeneratedParticles(n,m).resize(tmp);
        }
    }

    for(int n = 0; n < Nphases; ++n)
    for(int m = 0; m < Nphases; ++m)
    for(unsigned int i = 0; i < GeneratedParticles(n,m).size(); ++i)
    {
        inp >> GeneratedParticles(n,m)[i].x  
            >> GeneratedParticles(n,m)[i].y  
            >> GeneratedParticles(n,m)[i].z;
        inp >> GeneratedParticles(n,m)[i].Q1 
            >> GeneratedParticles(n,m)[i].Q2 
            >> GeneratedParticles(n,m)[i].Q3;
        inp >> GeneratedParticles(n,m)[i].radius;
    }
    inp.close();
}
} //namespace openphase
