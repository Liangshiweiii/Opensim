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

#ifndef FLOWSOLVERLBM_H
#define FLOWSOLVERLBM_H

//#include "Base/Includes.h"
#include "Base/Node.h"
#include "GrainInfo.h"
#include "FluidDynamics/D3Q27.h"
namespace openphase
{

class Settings;
class PhaseField;
class Composition;
class Velocities;
class D3Q27;
class Node;
class BoundaryConditions;
class GrainInfo;

class FlowSolverLBM : public OPObject                                           ///<  Calculation of fluid flow and advective solute transport
{
public:
    using OPObject::ReadInput;
    FlowSolverLBM(void){};                                                      ///<  Empty constructor
    FlowSolverLBM(const Settings& locSettings,
            const std::string InputFileName = std::string(),
            size_t boundary = 1);                                               ///<  Read Settings and InputFile with constructor
    void Initialize(const Settings& locSettings,
            const std::string InputFileName = std::string(),
            size_t boundary = 1);                                               ///<  Allocates memory and initializes global settings
    void ReadInput(const std::string InputFileName);                            ///<  Reads input parameters from a file
    static inline double psi(double rho, double rho0)                           ///< density potential function
    {
    #ifdef DEBUG
            if (rho <= 0.0)
                throw std::invalid_argument("Negative or zero density");
    #endif
        return 1.0 - exp(-rho/rho0);
    };

    dVector3 CalculateFluidMomentum(const PhaseField& Phase) const;             ///< Calculates the momentum of the entire fluid in the system

    double BounceBack(const int i, const int j, const int k, const int ii,
            const int jj, const int kk, const int n, PhaseField& Phase,
            const BoundaryConditions& BC, double& H);                           ///<  BounceBack at Solid Interfaces
    std::vector<double> CalculateFluidMass(void) const;                         ///<  Calculates mass of fluid components
    size_t CountObstacleNodes(void) const;                                      ///<  Counts the number of obstacle nodes
    size_t CountFluidNodes(void) const;                                         ///<  Counts the number of fluid nodes
    virtual void SetObstacleNodes(const PhaseField& Phase);                     ///<  Sets density and momentum in obstacles to be used when solid moves
    void ApplyForces(PhaseField& Phase, const Composition& Cx,
            const Velocities& Vel, const BoundaryConditions& BC);               ///<  Applies forces
    void ApplyForces(PhaseField& Phase, const Velocities& Vel,
            const BoundaryConditions& BC);                                      ///<  Applies force without composition
    void CalcateDensityAndMomentum(void);                                       ///<  Calculates Density and momentum from lbDensity
    void CalculateFluidVelocities(const PhaseField& Phase, Velocities& Vel);    ///<  Calculates Fluid velocities
    void CalculateForceBenzi(const int i, const int j, const int k,
            PhaseField& Phase);                                                 ///<  Force contribution according to Benzi
    void CalculateForceBuoyancy(const int i, const int j, const int k,
            PhaseField& Phase, const Composition& Cx);                          ///<  Force contribution of Bouyancy
    void CalculateForceDrag(const int i, const int j, const int k,
            PhaseField& Phase, const Velocities& Vel);                          ///<  Force contribution of Gravitation
    void CalculateForceGravitation(const int i, const int j, const int k);      ///<  Force contribution of Gravitation
    void Collision();                                                           ///<  Processes the collisions. Adjusts center of mass properties of colliding particles
    void DetectObstacles(const PhaseField& Phase);                              ///<  Detects Obstacles (Detects change of obstacles for advection)
    void DetectObstaclesSimple(const PhaseField& Phase);                        ///<  Detects Obstacles
    void EnforceSolidMomentum(PhaseField& Phase, const dVector3 value);         ///<  Enforces solid momentum conservation
    void EnforceMassConservation(void);                                         ///<  Enforces fluid mass conservation for advectional problems
    void FixPopulations(void);                                                  ///<  Ensures positive lbPopulations

    void Propagation(PhaseField& Phase, const BoundaryConditions& BC);          ///<  Propagates Populations
    void Read(const int tStep);                                                 ///<  Read raw (binary) fields from file
    void SetBoundaryConditions(const BoundaryConditions& BC);                   ///<  Sets boundary conditions for the particle distribution functions and flow velocities
    void SetFluidNodesNearObstacle();                                           ///<  Sets lbPoblulations of vanishing Obstacles and of appearing Obstacles
    void SetUniformVelocity(const BoundaryConditions& BC, const dVector3 U0);   ///<  Sets the initial values of particle distribution functions
    void Solve(PhaseField& Phase, const Composition& Cx, Velocities& Vel,
            const BoundaryConditions& BC);                                      ///<  Calculates one time step of the Navier-Stokes solver
    void Solve(PhaseField& Phase, Velocities& Vel,
            const BoundaryConditions& BC);                                      ///<  Calculates one time step of the Navier-Stokes solver
    void Write(const int tStep) const;                                          ///<  Write raw (binary) fields to file
    void WriteVTK(const int tStep) const;                                       ///<  Write fields into file with VTK format

    void Remesh(const int newNx, const int newNy, const int newNz,
            const BoundaryConditions& BC);

    Storage3D< D3Q27,    1 > lbPopulations;                                     ///<  Populations (discretized particle distribution functions) PDF)
    Storage3D< D3Q27,    1 > lbPopulationsTMP;                                  ///<  Temporary array for Populations PDF propagation
    Storage3D< bool,     0 > Obstacle;                                          ///<  1 if Node is solid
    Storage3D< bool,     0 > ObstacleAppeared;                                  ///<  True if obstacle note appeared
    Storage3D< bool,     0 > ObstacleVanished;                                  ///<  True if obstacle note vanished
    Storage3D< dVector3, 1 > lbForceDensity;                                    ///<  Force
    Storage3D< dVector3, 1 > lbMomentum;                                        ///<  Momentum
    Storage3D< dVector3, 1 > lbMomentumDelta;                                   ///<  Momentum change due to action of forces
    Storage3D< double,   1 > lbDensity;                                         ///<  Fluid density in lattice units

    bool ObstaclesChanged;                                                      ///<  True if obstacles changed

    dVector3 GA;                                                                ///<  Gravity acceleration in true units
    dVector3 lbGA;                                                              ///<  Gravity acceleration in lattice units
    double Eta;                                                                 ///<  Interface width in physical units
    //double RhoLiq;  //UNUSED                                                            ///<  Average liquid density
    //double RhoSol;  //UNUSED                                                            ///<  Average solid density
    double dRho;                                                                ///<  Coefficient to couple physical and lattice variables
    double dt;                                                                  ///<  Time step
    double dtx;                                                                 ///<  Coefficient to couple physical and lattice variables
    double dx3;                                                                 ///<  dx*dx*dx
    double dx;                                                                  ///<  Grid spacing
    double dxt2;                                                                ///<  dx/dt/dt
    double dxt;                                                                 ///<  Coefficient to couple lattice and physical variables
    double eta;                                                                 ///<  Physical interface width
    double h_star;                                                              ///<  Solid-liquid interaction parameter
    int N_Comp;                                                                 ///<  Number of components
    int N_Fluid_Comp;                                                           ///<  Number of fluid components
    int Nphases;                                                                ///<  Number of thermodynamic phases per grain
    int Nx;                                                                     ///<  Size of the inner calculation domain along X
    int Ny;                                                                     ///<  Size of the inner calculation domain along Y
    int Nz;                                                                     ///<  Size of the inner calculation domain along Z
    std::vector<double> FluidMass;                                              ///<  Fluid mass (used for mass conservation)
    //std::vector<double> V_star; //UNUSED                                                ///<  Solid-liquid interaction parameter
    std::vector<double> drhodc;                                                 ///<  Density change according to the change of concentration
    std::vector<double> nu;                                                     ///<  Kinematic viscosity
    std::vector<double> rho_0;                                                  ///<  Reference density (Parameter for Benzi)
    std::vector<double> tau;                                                    ///<  Relaxation time for BGK collision operator
    std::vector<std::vector<double> > G;                                        ///<  Parameter for Benzi
    long int FluidNodes;                                                        ///<  Stores number of obstacle nodes for next iteration

    bool Do_Benzi;                                                              ///<  Set to "true" if Benzi force should be calculated
    bool Do_BounceBack;                                                         ///<  Set to "true" if drag force at the interfaces should be calculated
    bool Do_Buoyancy;                                                           ///<  Set to "true" if buoyancy force should be calculated
    bool Do_Drag;                                                               ///<  Set to "true" if drag force at the interfaces should be calculated
    bool Do_Gravitation;                                                        ///<  Set to "true" if gravitational force should be calculated
    bool Do_SolidSolid;                                                         ///<  Set to "true" if drag force at the interfaces should be calculated
protected:
};

} //namespace openphase
#endif
