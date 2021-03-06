#include "Tools/Node.h"
#include "BoundaryConditions.h"
#include "GrainInfo.h"
#include "InterfaceMobility.h"
#include "PhaseField.h"
#include "PhysicalConstants.h"
#include "Temperatures.h"

namespace opensim
{
using namespace std;
/*************************************************************************/
InterfaceMobility::InterfaceMobility(const Settings& locSettings,
        unsigned int boundary)
{
    Initialize(locSettings, boundary);
    ReadInput(DefaultInputFileName);
}

void InterfaceMobility::Initialize(const Settings& locSettings, unsigned int boundary)
{
    thisclassname = "InterfaceMobility";
    //DefaultInputFileName = ProjectInputDir + "IntMobilityInput.opi";

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;
    Nphases = locSettings.Nphases;

    R = PhysicalConstants::R;

    Averaging = true;
    iWidth = locSettings.iWidth;

    IntMobility3D.Allocate(Nx, Ny, Nz, boundary);
    RawIntMobility3D.Allocate(Nx, Ny, Nz, boundary);
    Counter.Allocate(Nx, Ny, Nz, boundary);

    IntMobility.Allocate(Nphases, Nphases);
    IntAnisotropy.Allocate(Nphases, Nphases);
    MaxIntMobility.Allocate(Nphases, Nphases);
    ActivationEnergy.Allocate(Nphases, Nphases);

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void InterfaceMobility::ReadInput(const string InputFileName)
{
//  Read values of interface mobility for pairs of phases

    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname);
        exit(1);
    }
    Info::WriteBlankLine();
    Info::WriteLineInsert("Interface mobility");
    Info::WriteStandard("Source", InputFileName.c_str());
    
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
    
    for (int alpha = 0; alpha < Nphases; alpha++)
    for (int beta = alpha; beta < Nphases; beta++)
    {
        stringstream converter;
        converter << alpha << "_" << beta;
        string counter = converter.str();
        IntMobility(alpha, beta)      = UserInterface::ReadParameterD(inp, moduleLocation, string("Mu_") + counter);
        IntMobility(beta, alpha)      = IntMobility(alpha, beta);

        ActivationEnergy(alpha, beta) = UserInterface::ReadParameterD(inp, moduleLocation, string("AE_") + counter);
        ActivationEnergy(beta, alpha) = ActivationEnergy(alpha, beta);

        IntAnisotropy(alpha, beta)    = UserInterface::ReadParameterD(inp, moduleLocation, string("Eps_") + counter);
        IntAnisotropy(beta, alpha)    = IntAnisotropy(alpha, beta);
    }

    inp.close();

    GrainBoundariesAreMobile = false;

    for (int alpha = 0; alpha < Nphases; alpha++)
    if (IntMobility(alpha, alpha) > 0.0)
    {
        GrainBoundariesAreMobile = true;
        Info::WriteStandard(thisclassname, "Grain boundaries are mobile");
        break;
    }

    for(int n = 0; n < Nphases; n++)
    for(int m = 0; m < Nphases; m++)
    {
        MaxIntMobility(n,m) = IntMobility(n, m)*max(1.0, (1.0 + IntAnisotropy(n,m)));
    }
}

void InterfaceMobility::ReInitialize(const PhaseField& Phase)
{
    Nx = Phase.Nx;
    Ny = Phase.Ny;
    Nz = Phase.Nz;

    IntMobility3D.Reallocate(Nx, Ny, Nz);
    RawIntMobility3D.Reallocate(Nx, Ny, Nz);
    Counter.Reallocate(Nx, Ny, Nz);

    Info::WriteStandard(thisclassname, "Reinitialized");
}

void InterfaceMobility::Set(const PhaseField& Phase)
{
    double locmaxMu = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(max:locmaxMu))
    {
        if (Phase.Fields(i,j,k).flag)
        {
            clear(i,j,k);
            
            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
            for(auto  beta = alpha + 1;
                      beta != Phase.Fields(i,j,k).cend(); ++beta)
            {
                int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                double locIntMobility = IntMobility(pIndexA, pIndexB);
                set(i,j,k,alpha->index, beta->index, locIntMobility);
                locmaxMu = max(locIntMobility,locmaxMu);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxMu = locmaxMu;
}

void InterfaceMobility::Set(const PhaseField& Phase, const Temperature& Tx)
{
    double locmaxMu = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(max:locmaxMu))
    {
        if (Phase.Fields(i,j,k).flag)
        {
            clear(i,j,k);
            double invT = 1.0/(R*Tx(i,j,k));

            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
            for(auto  beta = alpha + 1;
                      beta != Phase.Fields(i,j,k).cend(); ++beta)
            {
                int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                double locIntMobility = IntMobility(pIndexA, pIndexB)*exp(-ActivationEnergy(pIndexA, pIndexB)*invT);
                set(i,j,k,alpha->index, beta->index, locIntMobility);
                locmaxMu = max(locIntMobility,locmaxMu);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxMu = locmaxMu;
}

void InterfaceMobility::CalculateCubic(const PhaseField& Phase)
{
    double locmaxMu = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(max:locmaxMu))
    {
        if (Phase.Fields(i,j,k).flag)
        {
            clear(i,j,k);

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

                    double MuCurr = IntMobility(pIndexA, pIndexB);
                    double AnisotropyCurr = IntAnisotropy(pIndexA, pIndexB);
                    double locMu = MuCurr*(1.0 - AnisotropyCurr * (1.5 - 2.5*
                                                         (pow(NormX, 4) +
                                                          pow(NormY, 4) +
                                                          pow(NormZ, 4))));

                    set(i,j,k,alpha->index, beta->index, locMu);
                    locmaxMu = max(locMu,locmaxMu);
                }
            }
            else
            {
                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta != Phase.Fields(i,j,k).cend(); ++beta)
                {
                    int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                    int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                    double locIntMobility = IntMobility(pIndexA, pIndexB);
                    set(i,j,k,alpha->index, beta->index, locIntMobility);
                    locmaxMu = max(locIntMobility,locmaxMu);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxMu = locmaxMu;
}

// void InterfaceMobility::CalculateHex(const PhaseField& Phase)
// {
//     double locmaxMu = 0.0;
//     OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(max:locmaxMu))
//     {
//         if (Phase.Fields(i,j,k).flag)
//         {
//             clear(i,j,k);

//             if (Phase.Interface(i,j,k))
//             {
//                 NodeV locNormals = Phase.Normals(i,j,k);

//                 for(auto alpha  = Phase.Fields(i,j,k).cbegin();
//                          alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
//                 for(auto  beta  = alpha + 1;
//                           beta != Phase.Fields(i,j,k).cend(); ++beta)
//                 {
//                     double NormX = 0.0;
//                     double NormY = 0.0;
//                     double NormZ = 0.0;

//                     locNormals.get(alpha->index, beta->index, NormX, NormY, NormZ);

//                     if(Phase.FieldsStatistics[alpha->index].State == Solid and 
//                        Phase.FieldsStatistics[ beta->index].State != Solid)
//                     {
//                         dVector3 Norm{NormX, NormY, NormZ};
//                         dVector3 NormR = Phase.FieldsStatistics[alpha->index].Orientation.RotationMatrix*Norm;

//                         NormX = NormR[0];
//                         NormY = NormR[1];
//                         NormZ = NormR[2];
//                     }
//                     if(Phase.FieldsStatistics[alpha->index].State != Solid and 
//                        Phase.FieldsStatistics[ beta->index].State == Solid)
//                     {
//                         dVector3 Norm{NormX, NormY, NormZ};
//                         dVector3 NormR = Phase.FieldsStatistics[ beta->index].Orientation.RotationMatrix*Norm;

//                         NormX = NormR[0];
//                         NormY = NormR[1];
//                         NormZ = NormR[2];
//                     }

//                     int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
//                     int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

//                     double MuCurr = IntMobility(pIndexA, pIndexB);
//                     double AnisotropyCurr = IntAnisotropy(pIndexA, pIndexB);

//                     double locIntMobility = MuCurr*(1.0 + AnisotropyCurr *
//                                          (      pow(NormX, 6) -
//                                                 pow(NormY, 6) -
//                                          15.0 * pow(NormX, 4) * NormY*NormY +
//                                          15.0 * pow(NormY, 4) * NormX*NormX +
//                                          (5.0 * pow(NormZ, 4) -
//                                           5.0 * pow(NormZ, 2) +
//                                                 pow(NormZ, 6))));

//                     set(i,j,k,alpha->index, beta->index, locIntMobility);
//                     locmaxMu = max(locIntMobility,locmaxMu);
//                 }
//             }
//             else
//             {
//                 for(auto alpha = Phase.Fields(i,j,k).cbegin();
//                          alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
//                 for(auto  beta = alpha + 1;
//                           beta != Phase.Fields(i,j,k).cend(); ++beta)
//                 {
//                     int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
//                     int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

//                     double locIntMobility = IntMobility(pIndexA, pIndexB);
//                     set(i,j,k,alpha->index, beta->index, locIntMobility);
//                     locmaxMu = max(locIntMobility,locmaxMu);
//                 }
//             }
//         }
//     }
//     OMP_PARALLEL_STORAGE_LOOP_END
//     maxMu = locmaxMu;
// }


void InterfaceMobility::CalculateCubic(const PhaseField& Phase, const Temperature& Tx)
{
    double locmaxMu = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(max:locmaxMu))
    {
        if (Phase.Fields(i,j,k).flag)
        {
            clear(i,j,k);

            if (Phase.Interface(i,j,k))
            {
                double invT = 1.0/(R*Tx(i,j,k));
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

                    double MuCurr = IntMobility(pIndexA, pIndexB)*exp(-ActivationEnergy(pIndexA, pIndexB)*invT);
                    double AnisotropyCurr = IntAnisotropy(pIndexA, pIndexB);
                    double locMu = MuCurr*(1.0 - AnisotropyCurr * (1.5 - 2.5*
                                                         (pow(NormX, 4) +
                                                          pow(NormY, 4) +
                                                          pow(NormZ, 4))));

                    set(i,j,k,alpha->index, beta->index, locMu);
                    locmaxMu = max(locMu,locmaxMu);
                }
            }
            else
            {
                double invT = 1.0/(R*Tx(i,j,k));

                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta != Phase.Fields(i,j,k).cend(); ++beta)
                {
                    int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                    int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                    double locIntMobility = IntMobility(pIndexA, pIndexB)*exp(-ActivationEnergy(pIndexA, pIndexB)*invT);
                    set(i,j,k,alpha->index, beta->index, locIntMobility);
                    locmaxMu = max(locIntMobility,locmaxMu);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxMu = locmaxMu;
}

void InterfaceMobility::CalculateHex(const PhaseField& Phase, const Temperature& Tx)
{
    double locmaxMu = 0.0;
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,reduction(max:locmaxMu))
    {
        if (Phase.Fields(i,j,k).flag)
        {
            clear(i,j,k);

            if (Phase.Interface(i,j,k))
            {
                double invT = 1.0/(R*Tx(i,j,k));
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

                    double MuCurr = IntMobility(pIndexA, pIndexB)*exp(-ActivationEnergy(pIndexA, pIndexB)*invT);
                    double AnisotropyCurr = IntAnisotropy(pIndexA, pIndexB);

                    double locIntMobility = MuCurr*(1.0 + AnisotropyCurr *
                                         (      pow(NormX, 6) -
                                                pow(NormY, 6) -
                                         15.0 * pow(NormX, 4) * NormY*NormY +
                                         15.0 * pow(NormY, 4) * NormX*NormX +
                                         (5.0 * pow(NormZ, 4) -
                                          5.0 * pow(NormZ, 2) +
                                                pow(NormZ, 6))));

                    set(i,j,k,alpha->index, beta->index, locIntMobility);
                    locmaxMu = max(locIntMobility,locmaxMu);
                }
            }
            else
            {
                double invT = 1.0/(R*Tx(i,j,k));

                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha != Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta != Phase.Fields(i,j,k).cend(); ++beta)
                {
                    int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                    int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                    double locIntMobility = IntMobility(pIndexA, pIndexB)*exp(-ActivationEnergy(pIndexA, pIndexB)*invT);
                    set(i,j,k,alpha->index, beta->index, locIntMobility);
                    locmaxMu = max(locIntMobility,locmaxMu);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxMu = locmaxMu;
}

void InterfaceMobility::Average(const PhaseField& Phase, const BoundaryConditions& BC)
{
    const int Range  = iWidth-1;
    const int Radius = iWidth-1;

    int tmpcellsPerLayer = 2*Range;
    int tmpslices = (Nx-1+tmpcellsPerLayer)/tmpcellsPerLayer;
    while (tmpslices % 2 == 1)
    {
        ++tmpcellsPerLayer;
        tmpslices = (Nx-1+tmpcellsPerLayer)/tmpcellsPerLayer;
    }
    const int cellsPerLayer = tmpcellsPerLayer;
    const int slices = tmpslices;
    const int Xrange = min(Range, Nx-1);
    const int Yrange = min(Range, Ny-1);
    const int Zrange = min(Range, Nz-1);
    if (Averaging)
    {
        #pragma omp parallel
        {
        int oddEven;
        #pragma omp master
        for (oddEven = 0; oddEven<=1; oddEven++)
        {
            int m;
            for (m=oddEven; m < slices; m+=2)
            #pragma omp task untied
            {
                Storage3D< int, 0 > SelectedPositions;
                SelectedPositions.Allocate(2*Xrange+1, 2*Yrange+1, 2*Zrange+1, 1);
                int max = std::min(cellsPerLayer * (m + 1),Nx);
                int i;
                for(i=cellsPerLayer * m; i < max; i++)
                {
                    for(int j = 0; j < Ny; j++)
                    for(int k = 0; k < Nz; k++)
                    if (Phase.Interface(i,j,k))
                    {
                        for(auto it = RawIntMobility3D(i,j,k).begin();
                                 it != RawIntMobility3D(i,j,k).end(); ++it)
                        {
                            double PhiAlpha = Phase.Fields(i,j,k).get(it->indexA);
                            double PhiBeta  = Phase.Fields(i,j,k).get(it->indexB);

                            if (PhiAlpha*PhiBeta > 0.09)
                            {
                                double DFSum = 0.0;
                                double SumWeights = 0.0;

                                for(int ii = -Xrange; ii <= Xrange; ii++)
                                for(int jj = -Yrange; jj <= Yrange; jj++)
                                for(int kk = -Zrange; kk <= Zrange; kk++)
                                /*if(i+ii >= 0 and i+ii < Nx and
                                   j+jj >= 0 and j+jj < Ny and
                                   k+kk >= 0 and k+kk < Nz)*/
                                {
                                    SelectedPositions(ii+Xrange, jj+Yrange, kk+Zrange) = 0.0;

                                    double dist = sqrt(ii*ii + jj*jj + kk*kk);

                                    double locPhiAlphaValue = Phase.Fields[BC.Index(i+ii, j+jj, k+kk, Nx, Ny, Nz, Phase.Fields.Bcells())][it->indexA];
                                    double  locPhiBetaValue = Phase.Fields[BC.Index(i+ii, j+jj, k+kk, Nx, Ny, Nz, Phase.Fields.Bcells())][it->indexB];

                                    double weight = (Radius - dist)*sqrt(locPhiAlphaValue * locPhiBetaValue);
                                    if (weight > DBL_EPSILON)
                                    {
                                        SumWeights += weight;
                                        DFSum += weight*RawIntMobility3D[BC.Index(i+ii, j+jj, k+kk, Nx, Ny, Nz, RawIntMobility3D.Bcells())].get_sym(it->indexA, it->indexB);
                                        SelectedPositions(ii+Xrange, jj+Yrange, kk+Zrange) = 1.0;
                                    }
                                }
                                double AveragedDF = 0.0;
                                if(SumWeights > 0.0) AveragedDF = DFSum/(SumWeights);

                                for(int ii = -Xrange; ii <= Xrange; ii++)
                                for(int jj = -Yrange; jj <= Yrange; jj++)
                                for(int kk = -Zrange; kk <= Zrange; kk++)
                                if(SelectedPositions(ii+Xrange, jj+Yrange, kk+Zrange) != 0.0)
                                {
                                    IntMobility3D[BC.Index(i+ii, j+jj, k+kk, Nx, Ny, Nz, IntMobility3D.Bcells())].add_sym(it->indexA, it->indexB, AveragedDF);
                                    Counter[BC.Index(i+ii, j+jj, k+kk, Nx, Ny, Nz, Counter.Bcells())].add_sym(it->indexA, it->indexB, 1.0);
                                }
                            }
                        }
                    }
                }
            #pragma omp taskwait
            } // odd/even
        }
        }
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
        if (Phase.Fields(i,j,k).flag)
        {
            if (IntMobility3D(i,j,k).size() > 0)
            {
                for(unsigned int n = 0; n < IntMobility3D(i,j,k).size(); ++n)
                {
                    auto averagedDF = (IntMobility3D(i,j,k).begin() + n);
                    auto counterDF = (Counter(i,j,k).begin() + n);
                    averagedDF->value /= counterDF->value;
                }
            }
            else
            {
                IntMobility3D(i,j,k) = RawIntMobility3D(i,j,k);
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    else
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
        if (Phase.Fields(i,j,k).flag)
        {
            IntMobility3D(i,j,k) = RawIntMobility3D(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

void InterfaceMobility::WriteVTK(const int tStep, const int indexA, const int indexB)
{
    stringstream converter;
    converter << indexA << "_" << indexB << "_" ;
    string phases = converter.str();
    string FileName = UserInterface::MakeFileName(VTKDir,"IntMobility_" + phases, tStep, ".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "IntMobility\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET RECTILINEAR_GRID\n";
    vtk_file << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    vtk_file << "X_COORDINATES " << Nx << " double\n";
    for (int i = 0; i < Nx; i++) vtk_file << i << " ";
    vtk_file << endl;
    vtk_file << "Y_COORDINATES " << Ny << " double\n";
    for (int i = 0; i < Ny; i++) vtk_file << i << " ";
    vtk_file << endl;
    vtk_file << "Z_COORDINATES " << Nz << " double\n";
    for (int i = 0; i < Nz; i++) vtk_file << i << " ";
    vtk_file << endl;
    vtk_file << "POINT_DATA " << Nx*Ny*Nz << "\n";

    vtk_file << "SCALARS Mu(" << indexA << "," << indexB << ") double" << "\n";
    vtk_file << "LOOKUP_TABLE default" << "\n";
    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        vtk_file << RawIntMobility3D(i,j,k).get_sym(indexA, indexB) << "\n";
    }

    vtk_file << "SCALARS MuAvg(" << indexA << "," << indexB << ") double" << "\n";
    vtk_file << "LOOKUP_TABLE default" << "\n";
    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        vtk_file << IntMobility3D(i,j,k).get_sym(indexA, indexB) << "\n";
    }
    vtk_file.close();
}

}// namespace opensim