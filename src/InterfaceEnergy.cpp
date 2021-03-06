#include "InterfaceEnergy.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "PhaseField.h"
#include "Tools/Node.h"
#include "GrainInfo.h"
#include "Tools.h"
#include "BoundaryConditions.h"

namespace opensim
{
using namespace std;
/*************************************************************************/
InterfaceEnergy::InterfaceEnergy(const Settings& locSettings,
        unsigned int boundary)
{
    Initialize(locSettings, boundary);
    ReadInput(DefaultInputFileName);
}

void InterfaceEnergy::Initialize(const Settings& locSettings,
        unsigned int boundary)
{
    thisclassname = "InterfaceEnergy";
    //DefaultInputFileName = ProjectInputDir + "IntEnergyInput.opi";

    boundary = max(boundary, (unsigned int)((locSettings.iWidth - 1)));//needed for driving force averaging

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;
    Nphases = locSettings.Nphases;
    IntEnergy3D.Allocate(Nx, Ny, Nz, boundary);
    RawIntEnergy3D.Allocate(Nx, Ny, Nz, boundary);
    Counter.Allocate(Nx, Ny, Nz, boundary);

    IntEnergy.Allocate(Nphases, Nphases);
    Anisotropy.Allocate(Nphases, Nphases);
    MinIntEnergy.Allocate(Nphases, Nphases);

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

void InterfaceEnergy::ReadInput(const string InputFileName)
{
    //  Read values of interface energy for pairs of phases

    fstream inp(InputFileName.c_str(), ios::in);
    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "Initialize");
        exit(1);
    }

    Info::WriteBlankLine();
    Info::WriteLineInsert("Interface energy anisotropic");
    Info::WriteStandard("Source", InputFileName);

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
    Info::WriteLine();
    inp.close();
}

void InterfaceEnergy::ReInitialize(const PhaseField& Phase)
{
    Nx = Phase.Nx;
    Ny = Phase.Ny;
    Nz = Phase.Nz;

    IntEnergy3D.Reallocate(Nx, Ny, Nz);

    Info::WriteStandard(thisclassname, "Reinitialized");
}

void InterfaceEnergy::Set(const PhaseField& Phase)
{
    double locmaxSigma = 0.0;
    const long int boundary = abs((Phase.Fields).Bcells()-1);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields, boundary,reduction(max:locmaxSigma))
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

                double locIntEnergy = IntEnergy(pIndexA, pIndexB);
                set(i,j,k,alpha->index, beta->index, locIntEnergy);
                locmaxSigma = max(locmaxSigma,locIntEnergy);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxSigma = locmaxSigma;
}

void InterfaceEnergy::SetRaw(const PhaseField& Phase)
{
    double locmaxSigma = 0.0;
    const long int boundary = abs((Phase.Fields).Bcells()-1);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields, boundary,reduction(max:locmaxSigma))
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

                double locIntEnergy = IntEnergy(pIndexA, pIndexB);
                set_raw(i,j,k,alpha->index, beta->index, locIntEnergy);
                locmaxSigma = max(locmaxSigma,locIntEnergy);
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxSigma = locmaxSigma;
}

void InterfaceEnergy::CalculateCubic(const PhaseField& Phase)
{
    double locmaxSigma = 0.0;
    const long int boundary = abs((Phase.Fields).Bcells()-1);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields, boundary,reduction(max:locmaxSigma))
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
                    double SigmaCurr = IntEnergy(pIndexA, pIndexB);
                    double AnisotropyCurr = Anisotropy(pIndexA, pIndexB);

                    double locIntEnergy = SigmaCurr*(1.0 + AnisotropyCurr *
                                            (1.5 - 2.5*(pow(NormX, 4) +
                                                        pow(NormY, 4) +
                                                        pow(NormZ, 4))));

                    set(i,j,k,alpha->index, beta->index, locIntEnergy);
                    locmaxSigma = max(locmaxSigma,locIntEnergy);
                }
            }
            else
            {
                for(auto alpha = Phase.Fields(i,j,k).cbegin();
                         alpha < Phase.Fields(i,j,k).cend() - 1; ++alpha)
                for(auto  beta = alpha + 1;
                          beta < Phase.Fields(i,j,k).cend(); ++beta)
                {
                    int pIndexA = Phase.FieldsStatistics[alpha->index].Phase;
                    int pIndexB = Phase.FieldsStatistics[ beta->index].Phase;

                    double locIntEnergy = MinIntEnergy(pIndexA, pIndexB);
                    set(i,j,k,alpha->index, beta->index, locIntEnergy);
                    locmaxSigma = max(locmaxSigma,locIntEnergy);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxSigma = locmaxSigma;
}

void InterfaceEnergy::CalculateHex(const PhaseField& Phase)
{
    double locmaxSigma = 0.0;
    const long int boundary = abs((Phase.Fields).Bcells()-1);
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,boundary,reduction(max:locmaxSigma))
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

                    set(i,j,k,alpha->index, beta->index, locIntEnergy);
                    locmaxSigma = max(locmaxSigma,locIntEnergy);
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

                    double locIntEnergy = MinIntEnergy(pIndexA, pIndexB);
                    set(i,j,k,alpha->index, beta->index, locIntEnergy);
                    locmaxSigma = max(locmaxSigma,locIntEnergy);
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxSigma = locmaxSigma;
}


void InterfaceEnergy::Average(const PhaseField& Phase, const BoundaryConditions& BC)
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
                        for(auto it = RawIntEnergy3D(i,j,k).begin();
                                 it != RawIntEnergy3D(i,j,k).end(); ++it)
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
                                        DFSum += weight*RawIntEnergy3D[BC.Index(i+ii, j+jj, k+kk, Nx, Ny, Nz, RawIntEnergy3D.Bcells())].get_sym(it->indexA, it->indexB);
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
                                	IntEnergy3D[BC.Index(i+ii, j+jj, k+kk, Nx, Ny, Nz, IntEnergy3D.Bcells())].add_sym(it->indexA, it->indexB, AveragedDF);
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
            if (IntEnergy3D(i,j,k).size() > 0)
            {
                for(unsigned int n = 0; n < IntEnergy3D(i,j,k).size(); ++n)
                {
                    auto averagedDF = (IntEnergy3D(i,j,k).begin() + n);
                    auto counterDF = (Counter(i,j,k).begin() + n);
                    averagedDF->value /= counterDF->value;
                }
            }
            else
            {
            	IntEnergy3D(i,j,k) = RawIntEnergy3D(i,j,k);
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    else
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
        if (Phase.Fields(i,j,k).flag)
        {
            IntEnergy3D(i,j,k) = RawIntEnergy3D(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
}

void InterfaceEnergy::WriteVTK(const int tStep, const int indexA,
        const int indexB)
{
    stringstream outbufer;

    outbufer << "# vtk DataFile Version 3.0\n";
    outbufer << "InterfaceEnergy\n";
    outbufer << "ASCII\n";
    outbufer << "DATASET STRUCTURED_GRID\n";
    outbufer << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    outbufer << "POINTS " <<  Nx*Ny*Nz << " double\n";

    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        outbufer << i << " " << j << " " << k << "\n";
    }
    outbufer << "\n";
    outbufer << "POINT_DATA " << Nx*Ny*Nz << "\n";

    outbufer << "SCALARS InterfaceEnergy_" << indexA
             << "_" << indexB << " double 1\n";
    outbufer << "LOOKUP_TABLE default\n";

    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        outbufer << IntEnergy3D(i,j,k).get_sym(indexA, indexB) << "\n";
    }

    stringstream tmp;
    tmp << "InterfaceEnergy_" << indexA << "_" << indexB << "_";
    string FileName = UserInterface::MakeFileName(VTKDir, tmp.str(), tStep, ".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbufer.rdbuf();
    vtk_file.close();
}  //  WriteVTK

void InterfaceEnergy::ConsiderStrain(const PhaseField& Phase,
        const ElasticProperties& EP)
{
    double locmaxSigma = 0.0;
    const long int boundary = abs((EP.EffectiveEigenStrains).Bcells());
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,EP.EffectiveEigenStrains,boundary,reduction(max:locmaxSigma))
    {
        const NodeV locNormals   = Phase.Normals(i,j,k);
        for (auto it = locNormals.cbegin(); it < locNormals.cend(); ++it)
        {
            // Calculate local projection matrix
            dMatrix3x3 Projection;
            Projection.set_to_zero();
            Projection(0,0) = 1 - it->X*it->X;
            Projection(1,0) =   - it->Y*it->X;
            Projection(2,0) =   - it->Z*it->X;
            Projection(0,1) =   - it->X*it->Y;
            Projection(1,1) = 1 - it->Y*it->Y;
            Projection(2,1) =   - it->Z*it->Y;
            Projection(0,2) =   - it->X*it->Z;
            Projection(1,2) =   - it->Y*it->Z;
            Projection(2,2) = 1 - it->Z*it->Z;

            double locCorrection = 1.0;
            for (int m = 0; m < 3; m++)
            for (int n = 0; n < 3; n++)
            {
                locCorrection += Projection(m,n) *
                    EP.Strains(i,j,k).get_tensor(m,n);
            }

            const double locIntEnergy =
                IntEnergy3D(i,j,k).get(it->indexA,it->indexB)*locCorrection;

            IntEnergy3D(i,j,k).set_sym(it->indexA,it->indexB,locIntEnergy);
            locmaxSigma = max(locmaxSigma,locIntEnergy);
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    maxSigma = locmaxSigma;
}

void InterfaceEnergy::WriteVTK(const int tStep,
        const PhaseField& Phase, const int indexA)
{
    stringstream outbufer;

    outbufer << "# vtk DataFile Version 3.0\n";
    outbufer << "InterfaceEnergy\n";
    outbufer << "ASCII\n";
    outbufer << "DATASET STRUCTURED_GRID\n";
    outbufer << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    outbufer << "POINTS " <<  Nx*Ny*Nz << " double\n";

    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        outbufer << i << " " << j << " " << k << "\n";
    }
    outbufer << "\n";
    outbufer << "POINT_DATA " << Nx*Ny*Nz << "\n";

    vector<int> presentPhaseFields;
    presentPhaseFields = Phase.GetPresentPhaseFields();

    for(unsigned int n = 0; n < presentPhaseFields.size(); n++)
    {
        int indexB = presentPhaseFields[n];

        bool phaseIsNeighbour = false;

        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            if(Phase.Fields(i,j,k)[indexA] > 0.0 and Phase.Fields(i,j,k)[indexB] > 0.0)
            {
                cout << "Phase " << indexA << "and phase " << indexB << " are neighbours:\n";
                Phase.PrintPointStatistics(i,j,k);
                phaseIsNeighbour = true;
            }
            if(phaseIsNeighbour)
                break;
        }

        if(phaseIsNeighbour)
        {
            outbufer << "SCALARS InterfaceEnergy_" << indexA
                     << "_" << indexB << " double 1\n";
            outbufer << "LOOKUP_TABLE default\n";

            for(int k = 0; k < Nz; ++k)
            for(int j = 0; j < Ny; ++j)
            for(int i = 0; i < Nx; ++i)
            {
                outbufer << IntEnergy3D(i,j,k).get_sym(indexA, indexB) << "\n";
            }
        }
    }

    stringstream tmp;
    tmp << "InterfaceEnergy_" << indexA << "_";
    string FileName = UserInterface::MakeFileName(VTKDir, tmp.str(), tStep, ".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbufer.rdbuf();
    vtk_file.close();
}  //  WriteVTK
}// namespace opensim