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

#ifndef NUCLEATION_H
#define NUCLEATION_H

#include "Tools/Includes.h"

namespace opensim
{
class PhaseField;
class Settings;
class InterfaceEnergy;
class DrivingForce;
class BoundaryConditions;
class Orientations;
class Temperature;
class SymmetryVariants;

class Nucleation : public OPObject                                           ///< Handles the nucleation of new phases and grains
{
  public:

    Nucleation(){};
    Nucleation(const Settings& locSettings, const std::string InputFileName = "NONE");

    void Initialize(const Settings& locSettings);                               ///< Initializes storages, sets internal variables.
    using OPObject::ReadInput;
    void ReadInput(const std::string InputFileName);                            ///< Reads input values from file
    void GenerateNucleationSites(PhaseField& Phase, Temperature& Tx);           ///< Randomly generates nucleation sites and gives them weights according to the chosen distribution parameters
    void ReGenerateNucleationSites(PhaseField& Phase, Temperature& Tx);         ///< Cleares existing seeds and generates new nucleation sites        
    void Clear();                                                               ///< Clears the particles storage
    void GenerateRandomSeeds(void);                                             ///< Generates new seeds for random number generators
    void PlantNuclei(PhaseField& Phi, int tstep);                               ///< Plants generated nuclei according to their nucleation parameters
    void PlantNuclei(PhaseField& Phi, SymmetryVariants& SV, int tstep);         ///< Plants generated nuclei according to their nucleation parameters
    void PlantNucleiGB(PhaseField& Phi, SymmetryVariants& SV, int tstep);       ///< Plants generated nuclei according to their nucleation parameters (GBs only)
    void CheckNuclei(PhaseField& Phi, InterfaceEnergy& IE, DrivingForce& dG, int tstep);///< Checks planted nuclei, removes unstable ones
    void CheckNuclei(PhaseField& Phi, InterfaceEnergy& IE, DrivingForce& dG, SymmetryVariants& SV, int tstep);///< Checks planted nuclei, removes unstable ones
    void WriteStatistics(int tstep, int PFindex, int NucleatingPhase, int MatrixPhase, 
                         int x, int y, int z, double Q1, double Q2, double Q3,
                         double dGnuc, double dGmin, std::string status) const; ///< Writes nulceation statistics to a file
    void Read(int tstep);                                                       ///< Read stored nucleation information from a file
    void Write(int tstep);                                                      ///< Write nucleation information to a file

    int  Nphases;                                                               ///< Number of thermodynamic phases
    int  Nx;                                                                    ///< X dimension of the system
    int  Ny;                                                                    ///< Y dimension of the system
    int  Nz;                                                                    ///< Z dimension of the system

    int SeedX;                                                                  ///< Seed for X dimension random number generator
    int SeedY;                                                                  ///< Seed for Y dimension random number generator
    int SeedZ;                                                                  ///< Seed for Z dimension random number generator
    
    // Parameters used to set particles size distribution for heterogeneous nucleation.
    // At the moment only normal distribution is implemented:
    // p = 1.0/(DistSigma*sqrt(2*Pi))*exp{-(x - DistMu)^2/(2.0*DistSigma^2)}
    
    struct Parameters                                                           ///< Pairwise nucleation settings ("phase alpha" in "phase beta")
    {
        double DistSigma;                                                       ///< Normal distribution standard deviation
        double DistMu;                                                          ///< Normal distribution mean value
        double Density;                                                         ///< Average nucleation density
        int    OrientationMode;                                                 ///< Nuclei orientation mode: Parent (inherites from the parent grain), Random (randomly generated) 
        int    Variants;                                                        ///< Nuclei symmetry variants selection: Single (randomly generated single variant), Multiple (generates multiple symmetry variants for each nucleation site)
        bool   Allowed;                                                         ///< True if nucleation is allowed in the bulk
        bool   AllowedInInterface;                                              ///< True if nucleation is allowed in the interface
        double Tmin;                                                            ///< Upper temperature limit for nucleation
        double Tmax;                                                            ///< Lower temperature limit for nucleation
        bool   Generated;                                                       ///< True if particles are already generated
        int    Mode;                                                            ///< Nucleation modes: bulk, grain boundary, both
    };

    Matrix< Parameters > NucleationParameters;                                  ///< Stores nucleation parameters for each pair of phases (A in B)
    
    struct NucSite                                                              ///< Particle properties
    {
        std::vector<int> PFindices;                                             ///< Phase-field indices after nucleation (for multiple symmetry variants)
        int PFindex;                                                            ///< Phase-field index after nucleation
        /// position indexes
        int x;
        int y;
        int z;
        /// orientation Euler angles
        double Q1;
        double Q2;
        double Q3;
        /// particle radius
        double radius;
        bool is_active()
        {
            return (not planted);
        }
        bool planted;                                                           ///< True if nucleus is planted successfully (and growing)
    };

    std::vector<std::vector<int> > GetNucleationSites(PhaseField& Phi);         ///< Vector of all nucleation events coordinates
    bool IsShielded(std::vector<NucSite> LocGrainStorage, const int i, const int j, const int k,
            const double shielding) const;                                      ///< Returns true if member of the "NucleatedGrains" is inside the "shielding"-radius, false otherwise.

    Matrix <std::vector <NucSite> > GeneratedParticles;                         ///< Generated particles storage
};

}// namespace opensim

#endif // NUCLEATION_H
