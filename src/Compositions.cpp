#include "Compositions.h"
#include "Info.h"
#include "Tools/Node.h"
#include "BoundaryConditions.h"
#include "Tools/UserInterface.h"
#include "VTK.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "Velocities.h"
#include "Chemistry/ChemicalProperties.h"
#include "Chemistry/PeriodicTable.h"

namespace opensim
{
    using namespace std;

    void Composition::Composition(const Settings& locSettings, const int boundary)
    {
        this->Initialize(locSettings, boundary);
        this->ReadInput(DefaultInputFileName);
    }

    void Composition::Initialize(const Settings& locSettings, const int boundarysize)
    {
        thisclassname = "Composition";
        //DefaultInputFileName = ProjectInputDir + "CompositionInput.opi";

        Nx      = locSettings.Nx;
        Ny      = locSettings.Ny;
        Nz      = locSettings.Nz;
        dx      = locSettings.dx;
        Ncomp   = locSettings.Ncomp-1;
        Nphases = locSettings.Nphases;
        AtStart = true;
        Threshold = DBL_EPSILON;
        Names.resize(Ncomp);
        MolarVolume.Allocate({Nphases, Ncomp});
        Min.Allocate({Nphases, Ncomp});
        Max.Allocate({Nphases, Ncomp});
        Phase.Allocate(Nx, Ny, Nz, {Nphases, Ncomp}, boundarysize);

        Norm.Allocate(Nx, Ny, Nz, {Nphases, Ncomp}, boundarysize);
        PhaseDot.Allocate(Nx, Ny, Nz, {Nphases, Ncomp}, boundarysize);

        Total.Allocate(Nx, Ny, Nz, {Ncomp}, boundarysize);
        Initial.Allocate({Nphases, Ncomp});
        Initial2.Allocate({Nphases, Ncomp});
        Limiting.Allocate(Nx, Ny, Nz, boundarysize);

        TotInitial.resize(Ncomp);

        initialized = true;
        Info::WriteStandard(thisclassname, "Initialized");
    }

    void Composition::Read(ChemicalProperties& CP, BoundaryConditions& BC,
                        int tStep, string File , bool legacy_format )
    {
        string FileName = UserInterface::MakeFileName(RawDataDir, File.c_str() , tStep, ".dat");

        fstream inp(FileName.c_str(), ios::in | ios::binary);

        if (!inp)
        {
            Info::WriteExit("File \"" + FileName + "\" could not be opened", thisclassname, "Read()");
            exit(1);
        };

        if(not legacy_format)
        {
            int locNx = Nx;
            int locNy = Ny;
            int locNz = Nz;
            int locNphases = Nphases;
            int locNcomp = Ncomp;
            inp.read(reinterpret_cast<char*>(&locNx), sizeof(int));
            inp.read(reinterpret_cast<char*>(&locNy), sizeof(int));
            inp.read(reinterpret_cast<char*>(&locNz), sizeof(int));
            inp.read(reinterpret_cast<char*>(&locNphases), sizeof(int));
            inp.read(reinterpret_cast<char*>(&locNcomp), sizeof(int));
            if(locNx != Nx or locNy != Ny or locNz != Nz or locNphases != Nphases or locNcomp != Ncomp)
            {
                stringstream message;
                message << "Inconsistent system dimensions!\n"
                        << "Input data dimensions: (" << locNx << ", " << locNy << ", " << locNz << ") grid points.\n"
                        << "Input data Nphases: " << locNphases << "\n"
                        << "Input data Ncomp: " << locNcomp << "\n"
                        << "Required data dimensions: (" << Nx << ", " << Ny << ", " << Nz << ") grid points.\n"
                        << "Required data Nphases: " << Nphases << "\n"
                        << "Required data Ncomp: " << Ncomp << "\n";
                Info::WriteExit(message.str(), thisclassname, "Read()");
                exit(1);
            }
        }
        STORAGE_LOOP_BEGIN(i,j,k,Phase,0)
            inp.read(reinterpret_cast<char*>(Phase(i,j,k).data()), Phase(i,j,k).size()*sizeof(double));
        STORAGE_LOOP_END
        STORAGE_LOOP_BEGIN(i,j,k,Total,0)
            inp.read(reinterpret_cast<char*>(Total(i,j,k).data()), Total(i,j,k).size()*sizeof(double));
        STORAGE_LOOP_END
        STORAGE_LOOP_BEGIN(i,j,k,Total,0)
            for(int n = 0; n < Ncomp; n++)
            {
                TotInitial[n] += Total(i, j, k)({n});
            }
        STORAGE_LOOP_END
        for(int n = 0; n < Ncomp; n++)
        {
            TotInitial[n] /= double(Nx*Ny*Nz);
        }
        SetBoundaryConditions(BC);
        CalculateTotalMolarVolume(CP);
        Info::WriteStandard(thisclassname, "Binary input loaded");
    }

    void Composition::WriteVTK(int tStep)
    {
        stringstream buffer;
        std::vector<int> DataTypes {PDScalars};

        VTK::WriteHeader(buffer, Nx, Ny, Nz);
        VTK::WriteBeginPointData(buffer, DataTypes);
        {
            for(int comp = 0; comp < Ncomp; comp++)
            {
                string nameComp = "TotalComposition_" + Names[comp];
                VTK::WriteVectorComponent(buffer, Total, nameComp, comp);
                for (int alpha = 0; alpha < Nphases; ++alpha)
                {
                    string namePhase = "PhaseComposition_" + Names[comp]
                                                    + "(" + to_string(alpha) + ")";
                    VTK::WriteTensor(buffer, Phase, namePhase, alpha, comp);
                }
            }
        }
        VTK::WriteEndPointData(buffer);
        VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
        VTK::WriteToFile(buffer, "Composition", tStep);
    }  //  WriteVTK

    void Composition::WriteVTK(PhaseField& Phi,ChemicalProperties& CP, Composition& Cx,int tStep)
    {
        /**This function will write VTK files for visualization of the 3D storages.
         * Total mole fractions, total weight percent, phase mole fractions as well
         * as site fractions.*/


        CalculateMoleFractions(CP);
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
            outbufer << "SCALARS TotalCompositionMF_" << CP.Component[comp].Name
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
            outbufer << "SCALARS TotalCompositionWP_" << CP.Component[comp].Name
                    << " double 1\n";
            outbufer << "LOOKUP_TABLE default\n";

            for(int x = 0; x < Nx; x++)
            for(int y = 0; y < Ny; y++)
            for(int z = 0; z < Nz; z++)
            {
                vector<double> WeightPercent(Ncomp, 0.0);
                WeightPercent = getWeightPercent(CP,Phi,x,y,z);

                outbufer << WeightPercent[comp] << "\n";
            }
        }

        for(int alpha = 0; alpha < Nphases; ++alpha)
        for(int comp = 0; comp < Ncomp; comp++)
        {
            outbufer << "SCALARS PhaseCompositionMF_" << CP.Component[comp].Name
                    << "(" << CP.Phase[alpha].Name << ") double 1\n";
            outbufer << "LOOKUP_TABLE default\n";

            for(int x = 0; x < Nx; x++)
            for(int y = 0; y < Ny; y++)
            for(int z = 0; z < Nz; z++)
            {
                if(Phi.Fractions(x,y,z)({alpha}) > 0.0)
                {
                    outbufer << get_MF(alpha,comp, x, y, z) << "\n";
                }
                else
                {
                    outbufer << 0.0 << "\n";
                }
            }
        }

        for(int alpha = 0; alpha < Nphases; ++alpha)
        {
            int index = 0;
            for(int sub = 0; sub < CP.Phase[alpha].Nsubs; ++sub)
            {
                for(int cons = 0; cons < CP.Phase[alpha].Sublattice[sub].Ncons; cons++)
                {
                    outbufer << "SCALARS PhaseCompositionSF_"
                            << CP.Phase[alpha].Sublattice[sub].Constituent[cons].Name
                            << "#" << sub << "(" << CP.Phase[alpha].Name
                            << ") double 1\n";
                    outbufer << "LOOKUP_TABLE default\n";

                    for(int x = 0; x < Nx; x++)
                    for(int y = 0; y < Ny; y++)
                    for(int z = 0; z < Nz; z++)
                    {
                        if(Phi.Fractions(x,y,z)({alpha}) > 0.0)
                        {
                            outbufer << get_CF(alpha,index,x,y,z) << "\n";
                        }
                        else
                        {
                            outbufer << 0.0 << "\n";
                        }
                    }

                    index++;
                }
            }
        }

        if(CP.Isotope)
        {
            for(int trac = 0; trac < CP.NumIsotopes; trac++)
            {
                outbufer << "SCALARS TotalTracerMF_" << CP.Component[CP.IsotopeDefined[trac]].Name
                        << " double 1\n";
                outbufer << "LOOKUP_TABLE default\n";

                for(int x = 0; x < Nx; x++)
                for(int y = 0; y < Ny; y++)
                for(int z = 0; z < Nz; z++)
                {
                    outbufer << CP.TotalTracer(x,y,z)({trac}) << "\n";
                }
            }
        }
        string FileName = UserInterface::MakeFileName(VTKDir, "Composition_",
                                                                tStep,".vtk");

        ofstream vtk_file(FileName.c_str());
        vtk_file << outbufer.rdbuf();
        vtk_file.close();

    }  //  WriteVTK

    void Composition::WriteStatistics(ChemicalProperties& CP, int tStep, double dt)
    {
        vector<double> total(CP.Ncomp, 0);
        vector<double> deviation(CP.Ncomp, 0);
        double sim_time = tStep*dt;

        for(int comp = 0; comp < CP.Ncomp; comp++)
        {
            for(int i = 0; i < CP.Nx; i++)
            for(int j = 0; j < CP.Ny; j++)
            for(int k = 0; k < CP.Nz; k++)
            {
                total[comp] += Total(i,j,k)({comp});
            }
            total[comp] /= double(CP.Nx*CP.Ny*CP.Nz);
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
            for(int comp = 0; comp < CP.Ncomp; comp++)
            {
                output_file << "total_" << CP.Component[comp].Name << "\t\t"
                            << "deviation_" << CP.Component[comp].Name << "\t\t";
            }
            output_file<< endl;
            output_file.close();
        }

        output_file.open("CompositionStatistics.txt", ios::app);
        output_file << tStep << "\t\t" << sim_time << "\t\t";
        for(int comp = 0; comp < CP.Ncomp; comp++)
        {
            output_file << total[comp]  << "\t\t" << deviation[comp] << "\t\t";
        }
        output_file << endl;
        output_file.close();
    }

    void Composition::WriteStatistics(int tStep, double dt)
    {
        vector<double> total(Ncomp, 0);
        vector<double> deviation(Ncomp, 0);
        double sim_time = tStep*dt;

        for(int comp = 0; comp < Ncomp; comp++)
        {
            for(int i = 0; i < Nx; i++)
            for(int j = 0; j < Ny; j++)
            for(int k = 0; k < Nz; k++)
            {
                total[comp] += Total(i,j,k)({comp});
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
                output_file << "total_" << Names[comp] << "\t\t"
                            << "deviation_" << Names[comp] << "\t\t";
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

    void Composition::Write(int tStep,  const string File, bool legacy_format)
    {
        string FileName = UserInterface::MakeFileName(RawDataDir, File.c_str(), tStep, ".dat");

        ofstream out(FileName.c_str(), ios::out | ios::binary);

        if (!out)
        {
            Info::WriteExit("File \"" + FileName + "\" could not be created", thisclassname, "Write()");
            exit(1);
        };

        if(not legacy_format)
        {
            int tmp = Nx;
            out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
            tmp = Ny;
            out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
            tmp = Nz;
            out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
            tmp = Nphases;
            out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
            tmp = Ncomp;
            out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
        }

        STORAGE_LOOP_BEGIN(i,j,k,Phase,0)
            out.write(reinterpret_cast<char*>(Phase(i,j,k).data()), Phase(i,j,k).size()*sizeof(double));
        STORAGE_LOOP_END
        STORAGE_LOOP_BEGIN(i,j,k,Total,0)
            out.write(reinterpret_cast<char*>(Total(i,j,k).data()), Total(i,j,k).size()*sizeof(double));
        STORAGE_LOOP_END
    }

    
}