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

#ifndef ADVECTIONHR_H
#define ADVECTIONHR_H

#include "Base/CommonFunctions.h"
#include "Base/Includes.h"
#include "Base/UserInterface.h"
#include "BoundaryConditions.h"
#include "Composition.h"
#include "FluidDynamics/FlowSolverLBM.h"
#include "Info.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "Orientations.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Temperature.h"
#include "Temperature.h"
#include "Velocities.h"

namespace openphase
{
    template <typename T>
    struct has_size
    {
        static const int value = 0;
    };

    template <>
    struct has_size<double>
    {
        static const int value = 1;
    };

    template <>
    struct has_size<dVector3>
    {
        static const int value = 3;
    };

    template <>
    struct has_size<Quaternion>
    {
        static const int value = 4;
    };

    template <>
    struct has_size<vStress>
    {
        static const int value = 6;
    };

    inline double Slope_Limiter_Minmod(const double a, const double b)
    {
        return 0.5*((a > 0) - (a < 0)+(b > 0) - (b < 0))*std::min(std::fabs(a),std::fabs(b));
    }

    inline double Slope_Limiter_MC(const double a, const double b)
    {
        return 0.5*((a > 0) - (a < 0)+(b > 0) - (b < 0))*std::min(0.5*std::fabs(a+b),std::min(2.0*std::fabs(a),2.0*std::fabs(b)));
    }

    inline double Slope_Limiter_Superbee(const double a, const double b)
    {
        return 0.5*((a > 0) - (a < 0)+(b > 0) - (b < 0))*std::max(std::min(2.0*std::fabs(a),std::fabs(b)),std::min(std::fabs(a),2.0*std::fabs(b)));
    }

    inline double Slope_Limiter_Upwind(const double a, const double b)
    {
        return 0.0;
    }

template<class T, long int rank>
class AdvectionMethod
{
 public:
    void Advection(Storage3D<T, rank>& Field,
                   const BoundaryConditions& BC, const Velocities& Vel, const int direction,
                   const double dt, const double dx, double (*Limiter)(const double,const double))
    {
        Info::WriteExit("Advection of template Storage3D<T, rank> not yet implemented.",
                        "AdvectionMethod", "Advection()");
        raise(SIGABRT);
    }
};

template<class T>
class AdvectionMethod<T,0>
{
 public:
    void Advection(Storage3D<T, 0>& Field,
                   const BoundaryConditions& BC, const Velocities& Vel, const int direction,
                   const double dt, const double dx, double (*Limiter)(const double,const double))
    {
        const int Nx = Field.sizeX();
        const int Ny = Field.sizeY();
        const int Nz = Field.sizeZ();
        Storage3D<T, 0> FieldBackup;
        Storage3D<T, 0> FieldUpdate;
        FieldBackup.Allocate(Nx, Ny, Nz, Field.Bcells());
        FieldUpdate.Allocate(Nx, Ny, Nz, Field.Bcells());
        FieldBackup = Field;
        int CFL = 0;
        double locdt = dt;
        std::vector<int> dir(3);
        dir[0] = (direction == 0)?1:0;
        dir[1] = (direction == 1)?1:0;
        dir[2] = (direction == 2)?1:0;
        const int offset = Field.Bcells()-2;
        int its = 1;
        do {
            CFL = 0;
            for (int it = 0; it < its; ++it)
            {
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,offset,reduction(+:CFL))
                for (int n = 0; n < has_size<T>::value; ++n)
                {
                    double q   = Field(i,j,k)[n];
                    double qp  = Field(i+dir[0],j+dir[1],k+dir[2])[n];
                    double qm  = Field(i-dir[0],j-dir[1],k-dir[2])[n];
                    double qpp = Field(i+2*dir[0],j+2*dir[1],k+2*dir[2])[n];
                    double qmm = Field(i-2*dir[0],j-2*dir[1],k-2*dir[2])[n];
                    double v = Vel.Average(i,j,k)[direction];
                    double vp = Vel.Average(i+dir[0],j+dir[1],k+dir[2])[direction];
                    double vm = Vel.Average(i-dir[0],j-dir[1],k-dir[2])[direction];
                    double L0 = 0.5*Limiter(q-qm,qp-q);
                    double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
                    double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
                    double ap = 0.5*(vp+v);
                    double am = 0.5*(vm+v);
                    double Fp = std::max(ap,0.)*(q+L0)
                              + std::min(ap,0.)*(qp-Lp1);
                    double Fm = std::max(am,0.)*(qm+Lm1)
                              + std::min(am,0.)*(q-L0);
                    if (std::fabs(ap)*locdt > 0.49*dx || std::fabs(am)*locdt > 0.49*dx)
                    {
                        CFL = 1;
                    }
                    FieldUpdate(i,j,k)[n] = Field(i,j,k)[n] - locdt/dx*(Fp-Fm);
                }
                OMP_PARALLEL_STORAGE_LOOP_END
                if (CFL == 0)
                {
                    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,offset,)
                    for (int n = 0; n < has_size<T>::value; ++n)
                    {
                        Field(i,j,k)[n] = FieldUpdate(i,j,k)[n];
                    }
                    OMP_PARALLEL_STORAGE_LOOP_END
                    BC.SetX(Field);
                    BC.SetY(Field);
                    BC.SetZ(Field);
                }
                else
                {
                    locdt /= 2.;
                    Field = FieldBackup;
                }
            }
            its *= 2;
        }
        while (CFL>0);
    }
};

template<long int rank>
class AdvectionMethod<double,rank>
{
 public:
    void Advection(Storage3D<double, rank>& Field,
                   const BoundaryConditions& BC, const Velocities& Vel, const int direction,
                   const double dt, const double dx, double (*Limiter)(const double,const double))
    {
        const unsigned int totsize = Field(0,0,0).size();

        Storage3D<double, rank> FieldBackup;
        Storage3D<double, rank> FieldUpdate;
        FieldBackup.Allocate(Field);
        FieldUpdate.Allocate(Field);
        /*OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,2,)
            for (int n = 0; n < totsize; ++n) {
                FieldBackup(i,j,k)[n] = Field(i,j,k)[n];
            }
        OMP_PARALLEL_STORAGE_LOOP_END */
        FieldBackup = Field;
        int CFL = 0;
        double locdt = dt;
        std::vector<int> dir(3);
        dir[0] = (direction == 0)?1:0;
        dir[1] = (direction == 1)?1:0;
        dir[2] = (direction == 2)?1:0;
        const int offset = Field.Bcells()-2;
        int its = 1;
        do {
            CFL = 0;
            for (int it = 0; it < its; ++it)
            {
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,offset,reduction(+:CFL))
                for (unsigned int n = 0; n < totsize; ++n)
                {
                    double q   = Field(i,j,k)[n];
                    double qp  = Field(i+dir[0],j+dir[1],k+dir[2])[n];
                    double qm  = Field(i-dir[0],j-dir[1],k-dir[2])[n];
                    double qpp = Field(i+2*dir[0],j+2*dir[1],k+2*dir[2])[n];
                    double qmm = Field(i-2*dir[0],j-2*dir[1],k-2*dir[2])[n];
                    double v = Vel.Average(i,j,k)[direction];
                    double vp = Vel.Average(i+dir[0],j+dir[1],k+dir[2])[direction];
                    double vm = Vel.Average(i-dir[0],j-dir[1],k-dir[2])[direction];
                    double L0 = 0.5*Limiter(q-qm,qp-q);
                    double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
                    double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
                    double ap = 0.5*(vp+v);
                    double am = 0.5*(vm+v);
                    double Fp = std::max(ap,0.)*(q+L0)
                              + std::min(ap,0.)*(qp-Lp1);
                    double Fm = std::max(am,0.)*(qm+Lm1)
                              + std::min(am,0.)*(q-L0);
                    if (std::fabs(ap)*locdt > 0.49*dx || std::fabs(am)*locdt > 0.49*dx)
                    {
                        CFL = 1;
                    }
                    FieldUpdate(i,j,k)[n] = Field(i,j,k)[n] - locdt/dx*(Fp-Fm);
                }
                OMP_PARALLEL_STORAGE_LOOP_END
                if (CFL == 0)
                {
                    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,offset,)
                    for (unsigned int n = 0; n < totsize; ++n)
                    {
                        Field(i,j,k)[n] = FieldUpdate(i,j,k)[n];
                    }
                    OMP_PARALLEL_STORAGE_LOOP_END
                    BC.SetX(Field);
                    BC.SetY(Field);
                    BC.SetZ(Field);
                }
                else
                {
                    locdt /= 2.;
                    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,offset,)
                        for (unsigned int n = 0; n < totsize; ++n) {
                            Field(i,j,k)[n] = FieldBackup(i,j,k)[n];
                        }
                    OMP_PARALLEL_STORAGE_LOOP_END
                }
            }
            its *= 2;
        }
        while (CFL>0);
    }
};

template<>
class AdvectionMethod<double,0>
{
 public:
    void Advection(Storage3D<double, 0>& Field,
                   const BoundaryConditions& BC, const Velocities& Vel, const int direction,
                   const double dt, const double dx, double (*Limiter)(const double,const double))
    {
        const int Nx = Field.sizeX();
        const int Ny = Field.sizeY();
        const int Nz = Field.sizeZ();
        Storage3D<double, 0> FieldBackup;
        Storage3D<double, 0> FieldUpdate;
        FieldBackup.Allocate(Nx, Ny, Nz, Field.Bcells());
        FieldUpdate.Allocate(Nx, Ny, Nz, Field.Bcells());
        FieldBackup = Field;
        int CFL = 0;
        double locdt = dt;
        std::vector<int> dir(3);
        dir[0] = (direction == 0)?1:0;
        dir[1] = (direction == 1)?1:0;
        dir[2] = (direction == 2)?1:0;
        const int offset = Field.Bcells()-2;
        int its = 1;
        do {
            CFL = 0;
            for (int it = 0; it < its; ++it)
            {
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,offset,reduction(+:CFL))
                    double q   = Field(i,j,k);
                    double qp  = Field(i+dir[0],j+dir[1],k+dir[2]);
                    double qm  = Field(i-dir[0],j-dir[1],k-dir[2]);
                    double qpp = Field(i+2*dir[0],j+2*dir[1],k+2*dir[2]);
                    double qmm = Field(i-2*dir[0],j-2*dir[1],k-2*dir[2]);
                    double v = Vel.Average(i,j,k)[direction];
                    double vp = Vel.Average(i+dir[0],j+dir[1],k+dir[2])[direction];
                    double vm = Vel.Average(i-dir[0],j-dir[1],k-dir[2])[direction];
                    double L0 = 0.5*Limiter(q-qm,qp-q);
                    double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
                    double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
                    double ap = 0.5*(vp+v);
                    double am = 0.5*(vm+v);
                    double Fp = std::max(ap,0.)*(q+L0)
                              + std::min(ap,0.)*(qp-Lp1);
                    double Fm = std::max(am,0.)*(qm+Lm1)
                              + std::min(am,0.)*(q-L0);
                    if (std::fabs(ap)*locdt > 0.49*dx || std::fabs(am)*locdt > 0.49*dx)
                    {
                        CFL = 1;
                    }
                    FieldUpdate(i,j,k) = Field(i,j,k) - locdt/dx*(Fp-Fm);
                OMP_PARALLEL_STORAGE_LOOP_END
                if (CFL == 0)
                {
                    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,offset,)
                        Field(i,j,k) = FieldUpdate(i,j,k);
                    OMP_PARALLEL_STORAGE_LOOP_END
                    BC.SetX(Field);
                    BC.SetY(Field);
                    BC.SetZ(Field);
                }
                else
                {
                    locdt /= 2.;
                    Field = FieldBackup;
                }
            }
            its *= 2;
        }
        while (CFL>0);
    }
};



class AdvectionHR : OPObject
{
 public:

    AdvectionHR(){};
    AdvectionHR(Settings& locSettings, std::string InputFileName = "NONE")
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

    void Initialize(const Settings& Settings)
    {
        thisclassname = "AdvectionHR";
        //DefaultInputFileName = ProjectInputDir + "AdvectionInput.opi";

        scheme = 0;                                                                 // Standard scheme is Upwind

        Info::WriteStandard(thisclassname, "Initialized");
    }

    using OPObject::ReadInput;

    void ReadInput(const std::string InputFileName)
    {
        Info::WriteBlankLine();
        Info::WriteLineInsert(thisclassname);
        Info::WriteStandard("Source", InputFileName);

        std::fstream inp(InputFileName.c_str(), std::ios::in);

        if (!inp)
        {
            Info::WriteExit("File " + InputFileName + " could not be opened",
                            thisclassname, "ReadInput()");
            exit(1);
        };
        int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
        schemeString = UserInterface::ReadParameterS(inp, moduleLocation, "scheme");

        std::vector<std::string> schemeTypes;
        schemeTypes.push_back("Upwind");
        schemeTypes.push_back("Minmod");
        schemeTypes.push_back("MC");
        schemeTypes.push_back("Superbee");

        bool schemeFound = false;
        for (unsigned int i = 0; i < schemeTypes.size(); i++)
        if(schemeString == schemeTypes[i])
        {
            schemeFound = true;
            scheme = i;
        }

        if(not schemeFound)
        {
            Info::WriteWarning("No advection scheme specified in "
                               + InputFileName + "! Choosing upwind!", thisclassname,
                               "ReadInput()");
        }

        inp.close();
        Info::WriteLine();
    }

    template<class T, long int rank>
    void Advect(Storage3D<T, rank>& Field, const Velocities& Vel,
                                              const BoundaryConditions& BC, const double dx, const double dt, const int tStep)
    {
        if (Field.Bcells() < 2)
        {
            Info::WriteExit("Number of Bcells for storage needs to be 2 or higher.",
                            "AdvectionHR", "Advection(Field)");
            raise(SIGABRT);
        }
        AdvectField(Field,Vel,BC,dx,dt,tStep);
    }

    void Advect(PhaseField& Phase, const Velocities& Vel,
                                              const BoundaryConditions& BC, const double dx, const double dt, const int tStep)
    {
        if (Phase.Fields.Bcells() < 2)
        {
            Info::WriteExit("Number of Bcells for storage Phase.Fields needs to be 2 or higher.",
                            "AdvectionHR", "Advection(PhaseField)");
            raise(SIGABRT);
        }
        AdvectPhaseField(Phase, Vel, BC, dx, dt, tStep);
    }

    void Advect(PhaseField& Phase, FlowSolverLBM& LBM, const Velocities& Vel,
                                              const BoundaryConditions& BC, const double dx, const double dt, const int tStep)
    {
        if (Phase.Fields.Bcells() < 2)
        {
            Info::WriteExit("Number of Bcells for storage Phase.Fields needs to be 2 or higher.",
                            "AdvectionHR", "Advection(PhaseField)");
            raise(SIGABRT);
        }
        LBM.DetectObstaclesSimple(Phase);
        AdvectPhaseField(Phase, Vel, BC, dx, dt, tStep);
        LBM.DetectObstacles(Phase);
        LBM.SetFluidNodesNearObstacle();
    }
    void Advect(Composition& Cx, const Velocities& Vel,
                                              const BoundaryConditions& BC, const double dx, const double dt, const int tStep)
    {
        if (Cx.Phase.Bcells() < 2)
        {
            Info::WriteExit("Number of Bcells for storage Cx.Phase needs to be 2 or higher.",
                            "AdvectionHR", "Advection(Composition)");
            raise(SIGABRT);
        }
        if (Cx.Total.Bcells() < 2)
        {
            Info::WriteExit("Number of Bcells for storage Cx.Total needs to be 2 or higher.",
                            "AdvectionHR", "Advection(Composition)");
            raise(SIGABRT);
        }
        Advect(Cx.Phase, Vel, BC, dx, dt, tStep);
        Advect(Cx.Total, Vel, BC, dx, dt, tStep);
    }
    void Advect(Temperature& Tx, const Velocities& Vel,
                                              const BoundaryConditions& BC, const double dx, const double dt, const int tStep)
    {
        if (Tx.Tx.Bcells() < 2)
        {
            Info::WriteExit("Number of Bcells for storage Tx.Tx needs to be 2 or higher.",
                            "AdvectionHR", "Advection(Temperature)");
            raise(SIGABRT);
        }
        Advect(Tx.Tx, Vel, BC, dx, dt,tStep);
    }
//    void Advect(ElasticProperties& EP, const Velocities& Vel,
//                                              const BoundaryConditions& BC, const double dx, const double dt, const int tStep)
//    {
//        if (EP.AccumulatedStresses.Bcells() < 2)
//        {
//            Info::WriteExit("Number of Bcells for storage EP.AccumulatedStresses needs to be 2 or higher.",
//                            "AdvectionHR", "Advection(ElasticProperties)");
//            raise(SIGABRT);
//        }
//        Advect(EP.AccumulatedStresses, Vel, BC, dx, dt, tStep);
//    }
    void Advect(Orientations& OR, const Velocities& Vel,
                                              const BoundaryConditions& BC, const double dx, const double dt, const int tStep)
    {
        // Number of Bcells of OR.Quaternions is allocated with 1 in Orientations class. Hence, reallocate.
        if (OR.Quaternions.Bcells() < 2)
        {
            OR.Quaternions.SetNumBcells(2);
        }
        Advect(OR.Quaternions, Vel, BC, dx, dt, tStep);

    }

 private:

    template<class T, long int rank>
    void AdvectField(Storage3D<T, rank>& Field, const Velocities& Vel,
                                              const BoundaryConditions& BC, const double dx, const double dt, const int tStep)
    {
        AdvectionMethod<T,rank> Adv;
        BC.SetX(Field);
        BC.SetY(Field);
        BC.SetZ(Field);
        double (*Limiter)(const double,const double);
        switch(scheme)
        {
            case 1: // Minmod
            {
                Limiter = &Slope_Limiter_Minmod;
                break;
            }
            case 2: // Monotonized Central
            {
                Limiter = &Slope_Limiter_MC;
                break;
            }
            case 3: // Superbee
            {
                Limiter = &Slope_Limiter_Superbee;
                break;
            }
            default: // Upwind
            {
                Limiter = &Slope_Limiter_Upwind;
                break;
            }
        }
        if (tStep % 2 == 0)
        {
            Adv.Advection(Field, BC, Vel, 0, dt, dx, Limiter);
            Adv.Advection(Field, BC, Vel, 1, dt, dx, Limiter);
            Adv.Advection(Field, BC, Vel, 2, dt, dx, Limiter);
        }
        else
        {
            Adv.Advection(Field, BC, Vel, 2, dt, dx, Limiter);
            Adv.Advection(Field, BC, Vel, 1, dt, dx, Limiter);
            Adv.Advection(Field, BC, Vel, 0, dt, dx, Limiter);
        }
    }

    void Advection(PhaseField& Phase,
                   Storage3D<Node, 0>& PhaseFieldUpdate,
                   Storage3D<Node, 0>& PhaseFieldBackup,
                   const BoundaryConditions& BC, const Velocities& Vel, const int direction,
                   const double dt, const double dx, double (*Limiter)(const double,const double))
    {
        PhaseFieldUpdate.Clear();
        PhaseFieldBackup = Phase.Fields;
        int CFL = 0;
        double locdt = dt;
        std::vector<int> dir(3);
        dir[0] = (direction == 0)?1:0;
        dir[1] = (direction == 1)?1:0;
        dir[2] = (direction == 2)?1:0;
        const int offset = Phase.Fields.Bcells()-2;
        int its = 1;
        do {
            CFL = 0;
            for (int it = 0; it < its; ++it)
            {
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,offset,reduction(+:CFL))
                if (Phase.Fields(i,j,k).flag)
                {
                    size_t solids = 0;
                    for(auto it = Phase.Fields(i,j,k).cbegin();
                            it != Phase.Fields(i,j,k).cend(); it++)
                    if (Phase.FieldsStatistics[it->index].State == Solid)
                    {
                        solids++;
                    }

                    //if (solids <= 1)
                    //{
                    //    int pInd = 0;
                    //    for (auto alpha = Phase.Fields(i,j,k).cbegin(); alpha < Phase.Fields(i,j,k).cend(); alpha++)
                    //    {
                    //        if(Phase.FieldsStatistics[alpha->index].State == Solid)
                    //        {
                    //            pInd = alpha->index;
                    //        }
                    //    }

                    //    pInd = Phase.FieldsStatistics[pInd].Phase;

                    //    for (auto alpha = Phase.Fields(i,j,k).begin(); alpha < Phase.Fields(i,j,k).end(); alpha++)
                    //    {
                    //        int ind = alpha->index;
                    //        double q   = Phase.Fields(i,j,k)[ind];
                    //        double qp  = Phase.Fields(i+dir[0],j+dir[1],k+dir[2])[ind];
                    //        double qm  = Phase.Fields(i-dir[0],j-dir[1],k-dir[2])[ind];
                    //        double qpp = Phase.Fields(i+2*dir[0],j+2*dir[1],k+2*dir[2])[ind];
                    //        double qmm = Phase.Fields(i-2*dir[0],j-2*dir[1],k-2*dir[2])[ind];
                    //        double v  = Vel.Phase(i,j,k)({pInd})[direction];
                    //        double vp = Vel.Phase(i+dir[0],j+dir[1],k+dir[2])({pInd})[direction];
                    //        double vm = Vel.Phase(i-dir[0],j-dir[1],k-dir[2])({pInd})[direction];
                    //        double L0  = 0.5*Limiter(q-qm,qp-q);
                    //        double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
                    //        double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
                    //        double ap = 0.5*(vp+v);
                    //        double am = 0.5*(vm+v);
                    //        double Fp = std::max(ap,0.)*(q+L0)
                    //                  + std::min(ap,0.)*(qp-Lp1);
                    //        double Fm = std::max(am,0.)*(qm+Lm1)
                    //                  + std::min(am,0.)*(q-L0);
                    //        if (std::fabs(ap)*locdt > 0.49*dx || std::fabs(am)*locdt > 0.49*dx)
                    //        {
                    //            CFL = 1;
                    //        }
                    //        PhaseFieldUpdate(i,j,k).set(ind, Phase.Fields(i,j,k)[ind] - locdt/dx*(Fp-Fm));
                    //    }
                    //}
                    //else
                    //{
                        for (auto alpha = Phase.Fields(i,j,k).begin();
                                alpha < Phase.Fields(i,j,k).end(); alpha++)
                        if (Phase.FieldsStatistics[alpha->index].IsMobile == 1 &&
                            Phase.FieldsStatistics[alpha->index].State == Solid )
                        {
                            const size_t Nx = Phase.Nx;
                            const size_t Ny = Phase.Ny;
                            const size_t Nz = Phase.Nz;

                            dMatrix3x3 W;
                            W(0,0) = 0.0;
                            W(1,1) = 0.0;
                            W(2,2) = 0.0;
                            W(0,1) = -Phase.FieldsStatistics[alpha->index].aVel[2];
                            W(0,2) =  Phase.FieldsStatistics[alpha->index].aVel[1];
                            W(1,2) = -Phase.FieldsStatistics[alpha->index].aVel[0];
                            W(1,0) =  Phase.FieldsStatistics[alpha->index].aVel[2];
                            W(2,0) = -Phase.FieldsStatistics[alpha->index].aVel[1];
                            W(2,1) =  Phase.FieldsStatistics[alpha->index].aVel[0];
                            const dVector3 pos  = {double(i),double(j),double(k)};
                            const dVector3 posp = {double(i+dir[0]),double(j+dir[1]),double(k+dir[2])};
                            const dVector3 posm = {double(i-dir[0]),double(j-dir[1]),double(k-dir[2])};
                            dVector3 distanceCM;
                            dVector3 distanceCMp;
                            dVector3 distanceCMm;
                            CommonFunctions::CalculateDistancePeriodic(pos, Phase.FieldsStatistics[alpha->index].Rcm,distanceCM,  Nx, Ny, Nz);
                            CommonFunctions::CalculateDistancePeriodic(posp,Phase.FieldsStatistics[alpha->index].Rcm,distanceCMp, Nx, Ny, Nz);
                            CommonFunctions::CalculateDistancePeriodic(posm,Phase.FieldsStatistics[alpha->index].Rcm,distanceCMm, Nx, Ny, Nz);
                            const dVector3 R  = distanceCM *dx;
                            const dVector3 Rp = distanceCMp*dx;
                            const dVector3 Rm = distanceCMm*dx;
                            const dVector3 Vcm = Phase.FieldsStatistics[alpha->index].Vcm;
                            const dVector3 Velocity  = Vcm + W*R;
                            const dVector3 Velocityp = Vcm + W*Rp;
                            const dVector3 Velocitym = Vcm + W*Rm;

                            const int ind = alpha->index;
                            const double q   = Phase.Fields(i,j,k)[ind];
                            const double qp  = Phase.Fields(i+dir[0],j+dir[1],k+dir[2])[ind];
                            const double qm  = Phase.Fields(i-dir[0],j-dir[1],k-dir[2])[ind];
                            const double qpp = Phase.Fields(i+2*dir[0],j+2*dir[1],k+2*dir[2])[ind];
                            const double qmm = Phase.Fields(i-2*dir[0],j-2*dir[1],k-2*dir[2])[ind];
                            const double v   = Velocity[direction];
                            const double vp  = Velocityp[direction];
                            const double vm  = Velocitym[direction];
                            const double L0  = 0.5*Limiter(q-qm,qp-q);
                            const double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
                            const double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
                            const double ap  = 0.5*(vp+v);
                            const double am  = 0.5*(vm+v);
                            const double Fp  = std::max(ap,0.)*(q+L0)
                                      + std::min(ap,0.)*(qp-Lp1);
                            const double Fm  = std::max(am,0.)*(qm+Lm1)
                                      + std::min(am,0.)*(q-L0);
                            if (std::fabs(ap)*locdt > 0.49*dx || std::fabs(am)*locdt > 0.49*dx)
                            {
                                CFL = 1;
                            }
                            PhaseFieldUpdate(i,j,k).set(ind, Phase.Fields(i,j,k)[ind] - locdt/dx*(Fp-Fm));
                        }
                        else
                        {
                            const int ind = alpha->index;
                            const double q   = Phase.Fields(i,j,k)[ind];
                            const double qp  = Phase.Fields(i+dir[0],j+dir[1],k+dir[2])[ind];
                            const double qm  = Phase.Fields(i-dir[0],j-dir[1],k-dir[2])[ind];
                            const double qpp = Phase.Fields(i+2*dir[0],j+2*dir[1],k+2*dir[2])[ind];
                            const double qmm = Phase.Fields(i-2*dir[0],j-2*dir[1],k-2*dir[2])[ind];
                            const double v  = Vel.Phase(i,j,k)({ind})[direction];
                            const double vp = Vel.Phase(i+dir[0],j+dir[1],k+dir[2])({ind})[direction];
                            const double vm = Vel.Phase(i-dir[0],j-dir[1],k-dir[2])({ind})[direction];
                            const double L0  = 0.5*Limiter(q-qm,qp-q);
                            const double Lp1 = 0.5*Limiter(qp-q,qpp-qp);
                            const double Lm1 = 0.5*Limiter(qm-qmm,q-qm);
                            const double ap  = 0.5*(vp+v);
                            const double am  = 0.5*(vm+v);
                            const double Fp  = std::max(ap,0.)*(q+L0)
                                      + std::min(ap,0.)*(qp-Lp1);
                            const double Fm  = std::max(am,0.)*(qm+Lm1)
                                      + std::min(am,0.)*(q-L0);
                            if (std::fabs(ap)*locdt > 0.49*dx || std::fabs(am)*locdt > 0.49*dx)
                            {
                                CFL = 1;
                            }
                            PhaseFieldUpdate(i,j,k).set(ind, Phase.Fields(i,j,k)[ind] - locdt/dx*(Fp-Fm));
                        }
                    //}
                }
                OMP_PARALLEL_STORAGE_LOOP_END
                if (CFL == 0)
                {
                    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,offset,)
                    if (Phase.Fields(i,j,k).flag)
                    {
                        Phase.Fields(i,j,k) = PhaseFieldUpdate(i,j,k);
                    }
                    else
                    {
                        Phase.Fields(i,j,k).clear();
                        Phase.Fields(i,j,k).add(Phase.Fields(i,j,k).front().index, 1.0);
                    }
                    OMP_PARALLEL_STORAGE_LOOP_END
                    BC.SetX(Phase.Fields);
                    BC.SetY(Phase.Fields);
                    BC.SetZ(Phase.Fields);
                }
                else
                {
                    locdt /= 2.;
                    Phase.Fields = PhaseFieldBackup;
                }
            }
            its *= 2;
        }
        while (CFL>0);
    }

    void AdvectPhaseField(PhaseField& Phase,
            const Velocities& Vel, const BoundaryConditions& BC,
            const double dx, const double dt, const int tStep)
    {
        const int Nx = Phase.Fields.sizeX();
        const int Ny = Phase.Fields.sizeY();
        const int Nz = Phase.Fields.sizeZ();
        Storage3D<Node, 0> PhaseFieldBackup;
        Storage3D<Node, 0> PhaseFieldUpdate;
        PhaseFieldBackup.Allocate(Nx, Ny, Nz, Phase.Fields.Bcells());
        PhaseFieldUpdate.Allocate(Nx, Ny, Nz, Phase.Fields.Bcells());
        BC.SetX(Phase.Fields);
        BC.SetY(Phase.Fields);
        BC.SetZ(Phase.Fields);
        double (*Limiter)(const double,const double);
        switch(scheme)
        {
            case 1: // Minmod
            {
                Limiter = &Slope_Limiter_Minmod;
                break;
            }
            case 2: // Monotonized Central
            {
                Limiter = &Slope_Limiter_MC;
                break;
            }
            case 3: // Superbee
            {
                Limiter = &Slope_Limiter_Superbee;
                break;
            }
            default: // Upwind
            {
                Limiter = &Slope_Limiter_Upwind;
                break;
            }
        }
        if (tStep % 2 == 0)
        {
            Advection(Phase, PhaseFieldUpdate, PhaseFieldBackup, BC, Vel, 0, dt, dx, Limiter);
            Advection(Phase, PhaseFieldUpdate, PhaseFieldBackup, BC, Vel, 1, dt, dx, Limiter);
            Advection(Phase, PhaseFieldUpdate, PhaseFieldBackup, BC, Vel, 2, dt, dx, Limiter);
        }
        else
        {
            Advection(Phase, PhaseFieldUpdate, PhaseFieldBackup, BC, Vel, 2, dt, dx, Limiter);
            Advection(Phase, PhaseFieldUpdate, PhaseFieldBackup, BC, Vel, 1, dt, dx, Limiter);
            Advection(Phase, PhaseFieldUpdate, PhaseFieldBackup, BC, Vel, 0, dt, dx, Limiter);
        }

        // Do special kind of finalization
        const int offset = Phase.Fields.Bcells();
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,offset,)
        {
            // restrict the phase field values to the interface [0,1]
            if(Phase.Fields(i,j,k).flag)
            {
                for(auto it = Phase.Fields(i,j,k).begin();
                            it != Phase.Fields(i,j,k).end();)
                {
                    int erase = 0;
                    if (it->value >= 1.0 - DBL_EPSILON) it->value = 1.0;
                    if (it->value <= 0.0 + DBL_EPSILON)
                    {
                        it->value = 0.0; erase = 1;
                    }

                    if (erase)
                    {
                        it = Phase.Fields(i,j,k).erase(it);
                    }
                    else
                    {
                        ++it;
                    }
                }

                double total = 0.0;
                for(auto it = Phase.Fields(i,j,k).begin();
                            it != Phase.Fields(i,j,k).end(); it++)
                {
                    total += it->value;
                }


                // Determine number solid and fluid phase-field and the
                // respective solid and fluid fractions
                size_t FluidPhases = 0;
                //size_t SolidPhases = 0;
                double SolidFraction = 0;
                //double FluidFraction = 0;
                for(auto it = Phase.Fields(i,j,k).cbegin();
                        it != Phase.Fields(i,j,k).cend(); it++)
                if (Phase.FieldsStatistics[it->index].State != Solid)
                {
                    //FluidFraction += it->value;
                    FluidPhases++;
                }
                else
                {
                    SolidFraction += it->value;
                    //SolidPhases++;
                }

                if (SolidFraction > 1)
                {
                    for(auto it = Phase.Fields(i,j,k).begin();
                            it != Phase.Fields(i,j,k).end();)
                    if (Phase.FieldsStatistics[it->index].State != Solid)
                    {
                        it = Phase.Fields(i,j,k).erase(it);
                    }
                    else
                    {
                        ++it;
                    }
                }
                else if (FluidPhases > 0)
                {
                        //if there are solid and fluid phases-fields present
                        // normalize only the fluid phases
                        const double norm = (1.0 - total)/FluidPhases;
                        for(auto it = Phase.Fields(i,j,k).begin();
                                    it != Phase.Fields(i,j,k).end(); it++)
                        if (Phase.FieldsStatistics[it->index].State != Solid)
                        {
                            it->value += norm;
                            it->value2 = 0.0;
                        }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        //Phase.Finalize(BC,false);
        Phase.Finalize(BC);
    }

    std::string schemeString;
    int scheme;
};

}// namespace openphase
#endif
