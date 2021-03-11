#include "PhaseField.h"
#include "Info.h"
#include "Mechanics/ElasticitySolvers/ElasticitySolverSpectral.h"
#include "BoundaryConditions.h"
#include "Tools/UserInterface.h"
#include "VTK.h"
#include "Velocities.h"
#include "Chemistry/ChemicalProperties.h"
#include <map>
#include "FluidDynamics/FlowSolverLBM.h"
namespace opensim
{
    using namespace std;

    PhaseField::PhaseField(const Settings& locSettings, unsigned int boundary)
    {
        this->Initialize(locSettings, boundary);
    }

    int PhaseField::Initialize(const Settings& locSettings, unsigned int boundary, int df)
    {
        thisclassname = "PhaseField";
        resolution = 1;
        NucleationPresent = false;

        Nx = locSettings.Nx;
        Ny = locSettings.Ny;
        Nz = locSettings.Nz;

        dx  = locSettings.dx;
        iWidth = locSettings.iWidth;
        Eta = locSettings.Eta;

        Nphases = locSettings.Nphases;

        if (df) boundary = max(int(boundary), int((iWidth - 1)));//needed for driving force averaging

        Fields.Allocate(Nx, Ny, Nz, boundary);
        FieldsDot.Allocate(Nx, Ny, Nz, boundary);
        Fractions.Allocate(Nx, Ny, Nz, { Nphases }, boundary);

        #ifdef _OPENMP
        omp_locks.Allocate(Nx, Ny, Nz, Fields.Bcells());
        init_omp_locks();
        #endif

        RefVolume = Pi;

        if(Nx < iWidth) RefVolume *= (Nx+1)/2.0;
        else RefVolume *= 1.1*iWidth;
        if(Ny < iWidth) RefVolume *= (Ny+1)/2.0;
        else RefVolume *= 1.1*iWidth;
        if(Nz < iWidth) RefVolume *= (Nz+1)/2.0;
        else RefVolume *= 1.1*iWidth;

        initialized = true;
        Info::WriteStandard(thisclassname, "Initialized");
        return boundary;
    }

    void PhaseField::Read(const BoundaryConditions& Bc, int tStep)
    {
        string FileName =
            UserInterface::MakeFileName(RawDataDir,"PF_", tStep, ".dat");

        Read(FileName, Bc);
        FieldsStatistics.Read(tStep);
        Finalize(Bc);
    }

    void PhaseField::Read(string FileName, const BoundaryConditions& Bc)
    {
        fstream inp(FileName.c_str(), ios::in | ios::binary);

        if (!inp)
        {
            Info::WriteExit(FileName + " could not be opened",
                    thisclassname, "Read()");
            exit(1);
        };

        int maxPFindex = 0;

        int locNx = Nx;
        int locNy = Ny;
        int locNz = Nz;
        inp.read(reinterpret_cast<char*>(&locNx), sizeof(int));
        inp.read(reinterpret_cast<char*>(&locNy), sizeof(int));
        inp.read(reinterpret_cast<char*>(&locNz), sizeof(int));
        if(locNx != Nx or locNy != Ny or locNz != Nz)
        {
            stringstream message;
            message << "Inconsistent system dimensions!\n"
                    << "Input data dimensions: (" << locNx
                    << ", " << locNy << ", " << locNz << ") grid points.\n"
                    << "Required data dimensions: (" << Nx
                    << ", " << Ny << ", " << Nz << ") grid points.\n";
            Info::WriteExit(message.str(), thisclassname, "Read()");
            exit(1);
        }

        STORAGE_LOOP_BEGIN(i,j,k,Fields,0)
        {
            Fields(i,j,k).clear();
            Fields(i,j,k).flag = 0;
            int num = 0;
            inp.read(reinterpret_cast<char*>(&num), sizeof(int));                   // Fields(i,j,k).size()
            for(int n = 0; n < num; n++)
            {
                int   idx = 0;
                double val = 0.0;
                inp.read(reinterpret_cast<char*>(&idx), sizeof(int));               // Fields(i,j,k)->index
                inp.read(reinterpret_cast<char*>(&val), sizeof(double));            // Fields(i,j,k)->value
                Fields(i,j,k).set(idx, val);

                if (idx > maxPFindex)
                {
                    maxPFindex = idx;
                }
            }
        }
        STORAGE_LOOP_END

        inp.close();

        FieldsStatistics.Allocate(maxPFindex + 1);

        Finalize(Bc);
        Info::WriteStandard(thisclassname, "Binary input loaded");
    }
    //将相场信息返回到一个整数存储
    int PhaseField::AddGrainInfo(int PhaseIndex)
    {
        int locIndex = FieldsStatistics.add_grain(PhaseIndex);

        return locIndex;
    }

    void PhaseField::NormalizeIncrements(const BoundaryConditions& BC, const double dt)

    {
        /** This function limits phase-field increments for all present phase-field
        pairs, so that the actual phase-field values are within their natural 
        limits of 0.0 and 1.0.*/
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
        if (Fields(i,j,k).flag)
        {
            for(auto alpha = Fields(i,j,k).cbegin();
                    alpha != Fields(i,j,k).cend(); ++alpha)
            {
                /* Set clearly corrupted sets to zero */
                double dPsiAlphaPos = 0.0;
                double dPsiAlphaNeg = 0.0;
                for(auto beta = Fields(i,j,k).cbegin();
                        beta != Fields(i,j,k).cend(); ++beta)
                if(alpha != beta)
                {
                    double locdPsi = FieldsDot(i,j,k).get_asym(alpha->index,
                                                            beta->index);
                    if(locdPsi > 0.0) dPsiAlphaPos += locdPsi;
                    if(locdPsi < 0.0) dPsiAlphaNeg += locdPsi;
                }
                if((alpha->value == 0.0 and dPsiAlphaNeg < 0.0 and dPsiAlphaPos == 0.0)
                or (alpha->value == 1.0 and dPsiAlphaPos > 0.0 and dPsiAlphaNeg == 0.0))
                {
                    for(auto beta = Fields(i,j,k).cbegin();
                            beta != Fields(i,j,k).cend(); ++beta)
                    if(alpha != beta)
                    {
                        FieldsDot(i,j,k).set_sym(alpha->index, beta->index, 0.0);
                    }
                }
            }
            for(auto alpha = FieldsDot(i,j,k).begin();
                    alpha != FieldsDot(i,j,k).end();)
            {
                /* Remove zero-sets from the storage */
                if(alpha->value == 0.0)
                {
                    alpha = FieldsDot(i,j,k).erase(alpha);
                }
                else
                {
                    ++alpha;
                }
            }
            /* Limit increments! This is done in a while loop, to acknowledge all
            existing pair-contributions.*/
            int number_of_iterations = 0;
            bool LimitingNeeded = true;
            while (LimitingNeeded)
            {
                number_of_iterations++;
                LimitingNeeded = false;
                Node locIncrements;
                for(auto it = FieldsDot(i,j,k).cbegin();
                        it < FieldsDot(i,j,k).cend(); ++it)
                {
                    /* Collect increments */
                    if(it->value < 0.0)
                    {
                        locIncrements.add2(it->indexA, it->value);
                        locIncrements.add(it->indexB, -it->value);
                    }
                    if(it->value > 0.0)
                    {
                        locIncrements.add(it->indexA, it->value);
                        locIncrements.add2(it->indexB, -it->value);
                    }
                }
                Node locLimits;
                for(auto alpha = Fields(i,j,k).cbegin();
                        alpha != Fields(i,j,k).cend(); ++alpha)
                {
                    /* Calculate limits */
                    double posIncrement = locIncrements.get(alpha->index);
                    double negIncrement = locIncrements.get2(alpha->index);
                    double newPFvalue = alpha->value+(posIncrement+negIncrement)*dt;
                    locLimits.set2(alpha->index,1.0);
                    locLimits.set(alpha->index,1.0);
                    if(newPFvalue < 0.0)
                    {
                        double tmpLim = locLimits.get2(alpha->index);
                        double tmpLim2 = min(tmpLim,-(alpha->value+posIncrement*dt)
                                                    /(negIncrement*dt));
                        locLimits.set2(alpha->index, tmpLim2);
                    }
                    if(newPFvalue > 1.0)
                    {
                        double tmpLim = locLimits.get(alpha->index);
                        double tmpLim2 = min(tmpLim,(1.0-(alpha->value+negIncrement
                                                    *dt))/(posIncrement*dt));
                        locLimits.set(alpha->index,tmpLim2);
                    }
                }
                for(auto it = FieldsDot(i,j,k).begin();
                        it < FieldsDot(i,j,k).end(); ++it)
                {
                    /* Limit increments */
                    if(it->value < 0.0)
                    {
                        double tmpLim= min(locLimits.get2(it->indexA),
                                        locLimits.get(it->indexB));
                        it->value *= tmpLim;
                        if (tmpLim < 1.0)
                        LimitingNeeded = true;
                    }
                    if(it->value > 0.0)
                    {
                        double tmpLim = min(locLimits.get(it->indexA),
                                            locLimits.get2(it->indexB));
                        it->value *= tmpLim;
                        if (tmpLim < 1.0)
                        LimitingNeeded = true;
                    }
                }
            }
            for(auto alpha = FieldsDot(i,j,k).begin();
                    alpha != FieldsDot(i,j,k).end(); )
            {
                /* cleaning up values which are too small */
                if(fabs(alpha->value) < DBL_EPSILON)
                {
                    alpha = FieldsDot(i,j,k).erase(alpha);
                }
                else
                {
                    ++alpha;
                }
            }
            /* Skip limiting loop, if no convergence after 24 iterations*/
            if (number_of_iterations > 24)
            LimitingNeeded = false;
            /* End plausibility check */
            Node locField = Fields(i,j,k);
            for(auto psi = FieldsDot(i,j,k).cbegin();
                    psi < FieldsDot(i,j,k).cend(); ++psi)
            if(psi->value != 0.0)
            {
                locField.add(psi->indexA,  psi->value * dt);
                locField.add(psi->indexB, -psi->value * dt);
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        SetIncrementsBoundaryConditions(BC);
        CalculateFractions();
        
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
        if(Fields(i,j,k).flag)
        {
            double zero = 1.0;
            
            for(int n = 0; n < Nphases; n++)
            {
                double fraction = Fractions(i,j,k)({n});
                double newfraction = fraction;
                
                zero -= fraction;
                
                for(auto it = FieldsDot(i,j,k).begin();
                it != FieldsDot(i,j,k).end(); ++it)
                {
                    if((FieldsStatistics[it->indexA].Phase == n)
                    and(FieldsStatistics[it->indexB].Phase != n))
                    {
                        newfraction += it->value*dt;
                    }
                    else if((FieldsStatistics[it->indexA].Phase != n)
                        and(FieldsStatistics[it->indexB].Phase == n))
                    {
                        newfraction -= it->value*dt;
                    }
                }
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    

    void PhaseField::MergeIncrements(const BoundaryConditions& BC,
                                            const double dt, const bool finalize)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
        {
            if(Fields(i,j,k).flag)
            {
                for(auto psi = FieldsDot(i,j,k).cbegin();
                        psi < FieldsDot(i,j,k).cend(); ++psi)
                if(psi->value != 0.0)
                {
                    Fields(i,j,k).add(psi->indexA,  psi->value * dt);
                    Fields(i,j,k).add(psi->indexB, -psi->value * dt);
                }

                FieldsDot(i,j,k).clear();
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
        Finalize(BC, finalize);
    }

    void PhaseField::WriteVTK(const int tStep, const bool CurvatureOutput) const
    {
        stringstream buffer;

        VTK::WriteHeader(buffer, Nx, Ny, Nz);
        std::vector<int> DataTypes {PDScalars};
        VTK::WriteBeginPointData(buffer, DataTypes);
        {
            WriteVTKData(buffer, CurvatureOutput);
        }
        VTK::WriteEndPointData(buffer);
        VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
        VTK::WriteToFile(buffer, "PhaseField", tStep);
    }

    void PhaseField::Write(int tStep)
    {
        string FileName =
            UserInterface::MakeFileName(RawDataDir,"PF_", tStep, ".dat");

        fstream out(FileName.c_str(), ios::out | ios::binary);

        if (!out)
        {
            Info::WriteExit("File \"" + FileName + "\" could not be opened",
                    thisclassname, "Write()");
            exit(1);
        };

        int tmp = Nx;
        out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
        tmp = Ny;
        out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
        tmp = Nz;
        out.write(reinterpret_cast<char*>(&tmp), sizeof(int));

        STORAGE_LOOP_BEGIN(i,j,k,Fields,0)
        {
            int tmp = Fields(i,j,k).size();
            out.write(reinterpret_cast<char*>(&tmp), sizeof(int));
            for(auto n = Fields(i,j,k).cbegin();
                    n < Fields(i,j,k).cend(); ++n)
            {
                int idx = n->index;
                out.write(reinterpret_cast<char*>(&idx), sizeof(int));
                double val = n->value;
                out.write(reinterpret_cast<char*>(&val), sizeof(double));
            }
        }
        STORAGE_LOOP_END

        out.close();

        FieldsStatistics.Write(tStep);
    }

    void PhaseField::PrintPFVolumes() const
    {
        for(unsigned int idx = 0; idx < FieldsStatistics.size(); idx++)
        {
            if(FieldsStatistics[idx].Volume > 0.0)
            {
                Info::WriteStandardNarrow("PF", std::to_string(idx));
                Info::WriteStandardNarrow("Variant",
                        std::to_string(FieldsStatistics[idx].Variant));
                Info::WriteStandardNarrow("Stage",
                        std::to_string(FieldsStatistics[idx].Stage));
                Info::WriteStandardNarrow("Volume",
                        std::to_string(FieldsStatistics[idx].Volume));
            }
        }
    }

}