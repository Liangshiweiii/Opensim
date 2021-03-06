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

#ifndef INTERFACEENERGYSOLIDLIQUID_H
#define INTERFACEENERGYSOLIDLIQUID_H

#include "Tools/Includes.h"
namespace opensim
{
class PhaseField;
class Composition;
class Orientations;
class InterfaceEnergy;
class Node;

class InterfaceEnergySolidLiquid : public OPObject                              ///< Provides interface energy models for solid-liquid interfaces
{
 public:
    InterfaceEnergySolidLiquid(void){};                                         ///< Default constructor
    InterfaceEnergySolidLiquid(Settings& locSettings, std::string InputFileName = "NONE")
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
    };
    void Initialize(Settings& locSettings);                                     ///< Initializes internal storagies and variables
    using OPObject::ReadInput;
    void ReadInput(std::string InputFileName);                                  ///< Reads input from file



    double Lorentz(double angle, double la, double lb)                          //Lorentz Function
    {
        double ret = la*la/(la*la + (angle-lb)*(angle-lb));
        return ret;
    }
    double LorentzPP(double angle, double la, double lb)                        //second derivative of Lorentz with respect to ANGLE
    {
        double ret = 2.0*la*la*(-la*la + 3.0*(angle-lb)*(angle-lb)) / pow((la*la + (angle-lb)*(angle-lb)),3);
        return ret;
    }
//
//    double Pochhammer(double argument, int n)
//    {
//        double ret = argument;
//        for (int i=1;i<n;i++)
//        {
//            ret *= argument+i;
//        }
//        return ret;
//    }
//
//    struct FacetVector
//    {   double x;
//        double y;
//        double z;
//        double la;
//        double lc;
//    };
//
//    std::vector <FacetVector> Facets;                         /// Facet Storage


    void CalculateCubic(PhaseField& Phase, InterfaceEnergy& IE);                ///< Calculates anisotropic interface energy for cubic symmetry grains
    void CalculateHexBoettger(PhaseField& Phase, InterfaceEnergy& IE);          ///< Calculates anisotropic interface energy for hexagonal symmetry grains according to B.Böttger et al. Acta Mat 2006, previously known as "CalculateHex"
    void CalculateHexSun(PhaseField& Phase, InterfaceEnergy& IE);               ///< Calculates anisotropic interface energy for hexagonal symmetry grains according to D.Y.Sun et al. PRB 2006
    void CalculateHexYang(PhaseField& Phase, InterfaceEnergy& IE);              ///< Calculates anisotropic interface energy for hexagonal symmetry grains according to M.Yang et al. Acta Mat 2015
    void CalculateLorentz(PhaseField& Phase, InterfaceEnergy& IE);              ///< Calculates anisotropic interface energy for arbitrarily chosen Lorentzian/Cauchy casp-functions for interface energy minima, in order to achieve faceted crystals. pseudo-crystals possible.
    int Nphases;                                                                ///< Number of thermodynamic phases

    Matrix< double > IntEnergy;                                                 ///< Interface energy parameters for solid-liquid phase pairs
    Matrix< double > MinIntEnergy;                                              ///< Minimum inteface energy for solid-liquid phase pairs 
    Matrix< double > Anisotropy;                                                ///< Interface energy anisotropy parameters for solid-liquid phase pairs
};

}// namespace opensim
#endif
