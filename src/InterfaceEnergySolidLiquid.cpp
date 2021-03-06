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

#include "InterfaceEnergySolidLiquid.h"
#include "InterfaceEnergy.h"
#include "PhaseField.h"
#include "Tools/Node.h"
#include "GrainInfo.h"
#include "Tools.h"

#include "Info.h"
#include "Settings.h"
#include "Orientations.h"
#include "Compositions.h"

namespace opensim
{
using namespace std;
/*************************************************************************/

void InterfaceEnergySolidLiquid::Initialize(Settings& locSettings)
{
    thisclassname = "InterfaceEnergySolidLiquid";
    //DefaultInputFileName = ProjectInputDir + "IntEnergyInput.opi";

    Nphases = locSettings.Nphases;
    
    IntEnergy.Allocate(Nphases, Nphases);
    Anisotropy.Allocate(Nphases, Nphases);
    MinIntEnergy.Allocate(Nphases, Nphases);

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void InterfaceEnergySolidLiquid::ReadInput(string InputFileName)
{
    Info::WriteLineInsert("InterfaceEnergySolidLiquid properties");
    Info::WriteStandard("Source", InputFileName);
        
    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput");
        exit(1);
    }
    
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
    
    for (int alpha = 0; alpha < Nphases; alpha++)
    for (int beta  = alpha; beta < Nphases; beta++)
    {
        stringstream converter;
        converter << alpha << "_" << beta;
        string counter = converter.str();
        IntEnergy(alpha, beta)     = UserInterface::ReadParameterD(inp, moduleLocation, string("Sigma_") + counter);
        IntEnergy(beta, alpha)     = IntEnergy(alpha, beta);
        Anisotropy(alpha, beta)    = UserInterface::ReadParameterD(inp, moduleLocation, string("Eps_") + counter);
        Anisotropy(beta, alpha)    = Anisotropy(alpha, beta);
    }

    for(int n = 0; n < Nphases; n++)
    for(int m = 0; m < Nphases; m++)
    {
        MinIntEnergy(n,m) = IntEnergy(n, m)*min(1.0, (1.0 - Anisotropy(n,m)));
    }

    inp.close();
}

void InterfaceEnergySolidLiquid::CalculateCubic(PhaseField& Phase, InterfaceEnergy& IE)
{
    double locmaxSigma = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(max:locmaxSigma))
    {
        if (Phase.Fields(i,j,k).flag)
        {
            IE.clear(i,j,k);
            
            if (Phase.Interface(i,j,k))
            {
                NodeV locNormals = Phase.Normals(i,j,k);
                
                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta != Phase.Fields(i,j,k).cend(); ++beta)
                {
                    double NormX = 0.0;
                    double NormY = 0.0;
                    double NormZ = 0.0;

                    locNormals.get(alpha->index, beta->index, NormX, NormY, NormZ);

                    if(Phase.FieldsStatistics[alpha->index].State == Solid and 
                       Phase.FieldsStatistics[ beta->index].State != Solid)
                    {
                        dVector3 Norm{NormX, NormY, NormZ};
                        dVector3 NormR = Phase.FieldsStatistics[alpha->index].Orientation.RotationMatrix*Norm;

                        NormX = NormR[0];
                        NormY = NormR[1];
                        NormZ = NormR[2];
                    }
                    if(Phase.FieldsStatistics[alpha->index].State != Solid and 
                       Phase.FieldsStatistics[ beta->index].State == Solid)
                    {
                        dVector3 Norm{NormX, NormY, NormZ};
                        dVector3 NormR = Phase.FieldsStatistics[ beta->index].Orientation.RotationMatrix*Norm;

                        NormX = NormR[0];
                        NormY = NormR[1];
                        NormZ = NormR[2];
                    }

                    int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                    int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;
                    double SigmaCurr = IntEnergy(pIndexA, pIndexB);
                    double AnisotropyCurr = Anisotropy(pIndexA, pIndexB);

                    double locIntEnergy = SigmaCurr*(1.0 + AnisotropyCurr *
                                            (1.5 - 2.5*(pow(NormX, 4) +
                                                        pow(NormY, 4) +
                                                        pow(NormZ, 4))));

                    IE.set(i,j,k,alpha->index, beta->index, locIntEnergy);
                    locmaxSigma = max(locmaxSigma,locIntEnergy);
                }
            }
            else
            {
                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha < Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta < Phase.Fields(i,j,k).cend(); ++beta)
                if(Phase.FieldsStatistics[alpha->index].State +
                   Phase.FieldsStatistics[ beta->index].State == SolidLiquid)
                {
                    int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                    int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                    double locIntEnergy = MinIntEnergy(pIndexA, pIndexB);
                    IE.set(i,j,k,alpha->index, beta->index, locIntEnergy);
                    locmaxSigma = max(locmaxSigma,locIntEnergy);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    IE.maxSigma = locmaxSigma;
}

void InterfaceEnergySolidLiquid::CalculateHexBoettger(PhaseField& Phase, InterfaceEnergy& IE)
{
    double locmaxSigma = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(max:locmaxSigma))
    {
        if (Phase.Fields(i,j,k).flag)
        {
            IE.clear(i,j,k);
            
            if (Phase.Interface(i,j,k))
            {
                NodeV locNormals = Phase.Normals(i,j,k);
                
                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta != Phase.Fields(i,j,k).cend(); ++beta)
                {
                    double NormX = 0.0;
                    double NormY = 0.0;
                    double NormZ = 0.0;

                    locNormals.get(alpha->index, beta->index, NormX, NormY, NormZ);

                    if(Phase.FieldsStatistics[alpha->index].State == Solid and 
                       Phase.FieldsStatistics[ beta->index].State != Solid)
                    {
                        dVector3 Norm{NormX, NormY, NormZ};
                        dVector3 NormR = Phase.FieldsStatistics[alpha->index].Orientation.RotationMatrix*Norm;

                        NormX = NormR[0];
                        NormY = NormR[1];
                        NormZ = NormR[2];
                    }
                    if(Phase.FieldsStatistics[alpha->index].State != Solid and 
                       Phase.FieldsStatistics[ beta->index].State == Solid)
                    {
                        dVector3 Norm{NormX, NormY, NormZ};
                        dVector3 NormR = Phase.FieldsStatistics[ beta->index].Orientation.RotationMatrix*Norm;

                        NormX = NormR[0];
                        NormY = NormR[1];
                        NormZ = NormR[2];
                    }    
                    
                    int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                    int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;
                    double SigmaCurr = IntEnergy(pIndexA, pIndexB);
                    double AnisotropyCurr = Anisotropy(pIndexA, pIndexB);

                    double locIntEnergy = SigmaCurr*(1.0 - AnisotropyCurr *
                                         (      pow(NormX, 6) -
                                                pow(NormY, 6) -
                                         15.0 * pow(NormX, 4) * NormY*NormY +
                                         15.0 * pow(NormY, 4) * NormX*NormX +
                                         (5.0 * pow(NormZ, 4) -
                                          5.0 * pow(NormZ, 2) +
                                                pow(NormZ, 6))));

                    IE.set(i,j,k,alpha->index, beta->index, locIntEnergy);
                    locmaxSigma = max(locmaxSigma,locIntEnergy);
                }
            }            
            else
            {
                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta != Phase.Fields(i,j,k).cend(); ++beta)
                if(Phase.FieldsStatistics[alpha->index].State +
                   Phase.FieldsStatistics[ beta->index].State == SolidLiquid)
                {
                    int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                    int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                    double locIntEnergy = MinIntEnergy(pIndexA, pIndexB);
                    IE.set(i,j,k,alpha->index, beta->index, locIntEnergy);
                    locmaxSigma = max(locmaxSigma,locIntEnergy);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    IE.maxSigma = locmaxSigma;
}


void InterfaceEnergySolidLiquid::CalculateHexSun(PhaseField& Phase, InterfaceEnergy& IE)
{
    double locmaxSigma = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(max:locmaxSigma))
    {
        if (Phase.Fields(i,j,k).flag)
        {
            IE.clear(i,j,k);

            if (Phase.Interface(i,j,k))
            {
                NodeV locNormals = Phase.Normals(i,j,k);

                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta != Phase.Fields(i,j,k).cend(); ++beta)
                {
                    double NormX = 0.0;
                    double NormY = 0.0;
                    double NormZ = 0.0;

                    locNormals.get(alpha->index, beta->index, NormX, NormY, NormZ);

                    if(Phase.FieldsStatistics[alpha->index].State == Solid and 
                       Phase.FieldsStatistics[ beta->index].State != Solid)
                    {
                        dVector3 Norm{NormX, NormY, NormZ};
                        dVector3 NormR = Phase.FieldsStatistics[alpha->index].Orientation.RotationMatrix*Norm;

                        NormX = NormR[0];
                        NormY = NormR[1];
                        NormZ = NormR[2];
                    }
                    if(Phase.FieldsStatistics[alpha->index].State != Solid and 
                       Phase.FieldsStatistics[ beta->index].State == Solid)
                    {
                        dVector3 Norm{NormX, NormY, NormZ};
                        dVector3 NormR = Phase.FieldsStatistics[ beta->index].Orientation.RotationMatrix*Norm;

                        NormX = NormR[0];
                        NormY = NormR[1];
                        NormZ = NormR[2];
                    }            

                    int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                    int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;
                    double SigmaCurr = IntEnergy(pIndexA, pIndexB);
                    double eps20 = -0.026; //TO DO: Create Read In for HexSun
                    double eps40 =  0.0;
                    double eps60 =  0.0;
                    double eps66 =  0.003;

                    double locIntEnergy = SigmaCurr*(1.0
                            +eps20*sqrt(5.0/16.0/Pi)*(3.0*NormZ*NormZ - 1.0)
                            +eps40*3.0/16.0/sqrt(Pi)*(35.0*NormZ*NormZ*NormZ*NormZ - 30.0*NormZ*NormZ + 3.0)
                            +eps60*sqrt(13.0/Pi)/32.0*(231.0*NormZ*NormZ*NormZ*NormZ*NormZ*NormZ - 315.0*NormZ*NormZ*NormZ*NormZ + 105.0*NormZ*NormZ - 5.0)
                            +eps66*sqrt(6006.0/Pi)/64.0*(NormX*NormX*NormX*NormX*NormX*NormX - 15.0*NormX*NormX*NormX*NormX*NormY*NormY + 15.0*NormX*NormX*NormY*NormY*NormY*NormY - NormY*NormY*NormY*NormY*NormY*NormY));

                    IE.set(i,j,k,alpha->index, beta->index, locIntEnergy);
                    locmaxSigma = max(locmaxSigma,locIntEnergy);
                }
            }
            else
            {
                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta != Phase.Fields(i,j,k).cend(); ++beta)
                if(Phase.FieldsStatistics[alpha->index].State +
                   Phase.FieldsStatistics[ beta->index].State == SolidLiquid)
                {
                    int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                    int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                    double locIntEnergy = MinIntEnergy(pIndexA, pIndexB);
                    IE.set(i,j,k,alpha->index, beta->index, locIntEnergy);
                    locmaxSigma = max(locmaxSigma,locIntEnergy);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    IE.maxSigma = locmaxSigma;
}

void InterfaceEnergySolidLiquid::CalculateHexYang(PhaseField& Phase, InterfaceEnergy& IE)
{
    double locmaxSigma = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(max:locmaxSigma))
    {
        if (Phase.Fields(i,j,k).flag)
        {
            IE.clear(i,j,k);

            if (Phase.Interface(i,j,k))
            {
                NodeV locNormals = Phase.Normals(i,j,k);

                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta != Phase.Fields(i,j,k).cend(); ++beta)
                {
                    double NormX = 0.0;
                    double NormY = 0.0;
                    double NormZ = 0.0;

                    locNormals.get(alpha->index, beta->index, NormX, NormY, NormZ);

                    if(Phase.FieldsStatistics[alpha->index].State == Solid and 
                       Phase.FieldsStatistics[ beta->index].State != Solid)
                    {
                        dVector3 Norm{NormX, NormY, NormZ};
                        dVector3 NormR = Phase.FieldsStatistics[alpha->index].Orientation.RotationMatrix*Norm;

                        NormX = NormR[0];
                        NormY = NormR[1];
                        NormZ = NormR[2];
                    }
                    if(Phase.FieldsStatistics[alpha->index].State != Solid and 
                       Phase.FieldsStatistics[ beta->index].State == Solid)
                    {
                        dVector3 Norm{NormX, NormY, NormZ};
                        dVector3 NormR = Phase.FieldsStatistics[ beta->index].Orientation.RotationMatrix*Norm;

                        NormX = NormR[0];
                        NormY = NormR[1];
                        NormZ = NormR[2];
                    }
                    
                    int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                    int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;
                    double SigmaCurr = IntEnergy(pIndexA, pIndexB);
                    double eps1 = -0.02; //TO DO: Create Read In for HexYang
                    double eps2 =  0.15;
                    double eps3 =  0.15;

                    double locIntEnergy = SigmaCurr*(1.0
                            + eps1 * (3.0*NormZ*NormZ - 1.0)*(3.0*NormZ*NormZ - 1.0)
                            + eps2 * pow((NormX*NormX*NormX - 3.0*NormX*NormY*NormY),2)
                            * pow((9.0*NormZ*NormZ - 1.0 + eps3),2));

                    IE.set(i,j,k,alpha->index, beta->index, locIntEnergy);
                    locmaxSigma = max(locmaxSigma,locIntEnergy);
                }
            }
            else
            {
                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta != Phase.Fields(i,j,k).cend(); ++beta)
                if(Phase.FieldsStatistics[alpha->index].State +
                   Phase.FieldsStatistics[ beta->index].State == SolidLiquid)
                {
                    int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                    int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                    double locIntEnergy = MinIntEnergy(pIndexA, pIndexB);
                    IE.set(i,j,k,alpha->index, beta->index, locIntEnergy);
                    locmaxSigma = max(locmaxSigma,locIntEnergy);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    IE.maxSigma = locmaxSigma;
}

void InterfaceEnergySolidLiquid::CalculateLorentz(PhaseField& Phase, InterfaceEnergy& IE)
{
        //Using the grain orientation matrix, rotate the interface normal such that they are relative to that.
        //the result is the same, as if the grain was alligned with the parent coordinate system.
        //Interface-normal now consists of NormX, NormY, NormZ, RELATIVE to grain orientation.
        //NormZ equals cos(THETAone)
        //double cosThetaOneSqare = pow(NormZ,2);
        //double cosThetaTwo = NormX/sqrt(pow(NormX,2)+pow(NormY,2)); //Gx*Norm{notZ}/[abs(Gx)*abs(Norm{nozZ})] = cos(angPhi) , where angPhi is the angle between the vector Gx and Norm{notZ}. In our case Gx is set to (1,0,0) and Norm{notZ} := (NormX,NormY,0)
        //double ThetaTwo = acos(NormX/sqrt(pow(NormX,2)+pow(NormY,2)))
        //Lxy(ThetaTwo) := Lorentzian which will be used in the xy plane. has its peak at 0.
        //Lz(ThetaOne) := Lorentzian which will be used in the z related contribution. has its peak at 0.
//
//    for(int n = 0; n < Phase.Nphases; n++)
//    for(int m = 0; m < Phase.Nphases; m++)
//    {
//        MinIntEnergy(n,m) = IntEnergy(n, m)*0.05;
//    }

    double locmaxSigma = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(max:locmaxSigma))
    {
        if (Phase.Fields(i,j,k).flag)
        {
            IE.clear(i,j,k);

            if (Phase.Interface(i,j,k))
            {
                NodeV locNormals = Phase.Normals(i,j,k);

                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta != Phase.Fields(i,j,k).cend(); ++beta)
                {
                    double NormX = 0.0;
                    double NormY = 0.0;
                    double NormZ = 0.0;

                    locNormals.get(alpha->index, beta->index, NormX, NormY, NormZ);

                    if(Phase.FieldsStatistics[alpha->index].State == Solid and 
                       Phase.FieldsStatistics[ beta->index].State != Solid)
                    {
                        dVector3 Norm{NormX, NormY, NormZ};
                        dVector3 NormR = Phase.FieldsStatistics[alpha->index].Orientation.RotationMatrix*Norm;

                        NormX = NormR[0];
                        NormY = NormR[1];
                        NormZ = NormR[2];
                    }
                    if(Phase.FieldsStatistics[alpha->index].State != Solid and 
                       Phase.FieldsStatistics[ beta->index].State == Solid)
                    {
                        dVector3 Norm{NormX, NormY, NormZ};
                        dVector3 NormR = Phase.FieldsStatistics[ beta->index].Orientation.RotationMatrix*Norm;

                        NormX = NormR[0];
                        NormY = NormR[1];
                        NormZ = NormR[2];
                    }

//                for (auto it = Facets(x,y,z,la,lc).begin(); it < Facets(x,y,z,la,lc).end(); ++it)
//                {
//                    FacetVector locFacet;
//                    locFacet.x = XNorm;
//                    locFacet.y = YNorm;
//                    locFacet.z = Znorm;
//                    locFacet.la = LorentzA;
//                    locFacet.lc = Amplitude;
//
//                    Facets.push_back(locFacet);
//                }

                double ThetaOne = 0.0;
                double ThetaTwo = 0.0;
                ThetaOne = acos(NormZ)*180.0/Pi; //this is the angle between the Z axis and the normal. "theta"
                //cout <<"ThetaOne: " << ThetaOne << endl;
//                if (NormY >= 0.0)
//                    ThetaTwo = acos(NormX/sqrt(pow(NormX,2)+pow(NormY,2)))*180.0/Pi;
                //this is the angle between the X axis and the XY-plane normal-projection. "phi"
//                else
//                    ThetaTwo = 360.0 - acos(NormX/sqrt(pow(NormX,2)+pow(NormY,2)))*180.0/Pi;
                //cout <<"ThetaTwo: " << ThetaTwo << endl;
                ThetaTwo = atan2(NormY,NormX)*180.0/Pi+180.0;
                double la,lb; //parameters for Lorentzian "Lorentz(angle, a, b, c)"
                la = 5.1;
                lb = 0.0;
                //lc = 1.0;

                double casps   = 1.0;   //casps is intended to represent the product of all casps (Product over all (1-Lorentzians)) at the angle of the interface normal. EDIT: Now it is 1 minus eventually applicable casps.
                double caspsPP = 0.0; //caspsPP is intended to represent the sum of all casps of the second derivative of the Lorentzian at the angle of the interface normal.
                double ThetaOneSymmetry = 90.0;                                 //i.e. every ThetaOneSymmetry degrees, there is a casp.
                double ThetaTwoSymmetry = 72.0;


                for (int p=0; p<=180; p+=ThetaOneSymmetry)
                {
                    casps   -= Lorentz  (ThetaOne, la, lb+p);
                    if (casps <= 0.01) casps = 0.01;
                }

                for (double p=0; p<=360; p+=ThetaTwoSymmetry)
                {
                    casps   -= Lorentz  (ThetaTwo, la, lb+p);
                    if (casps <= 0.01) casps = 0.01;
                }

                for (double p=0; p<=180; p+=ThetaOneSymmetry)
                {
                    caspsPP += LorentzPP(ThetaOne, la, lb+p);
                }

                for (double p=0; p<=360; p+=ThetaTwoSymmetry)
                {
                    caspsPP += LorentzPP(ThetaTwo, la, lb+p);
                }

                int pIndexA = alpha->indexA % Nphases;
                int pIndexB = alpha->indexB % Nphases;
                double locIntEnergy = IntEnergy(pIndexA, pIndexB) * (casps+caspsPP);
                //cout << locIntEnergy << endl;
                IE.set(i,j,k,alpha->indexA, alpha->indexB, locIntEnergy);
                locmaxSigma = max(locmaxSigma,locIntEnergy);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    IE.maxSigma = locmaxSigma;
}
}// namespace openphase

