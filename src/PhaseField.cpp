
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

void PhaseField::Clear()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
    {
        Fields(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::CalculateFractions()
{
	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, Fields.Bcells(), )
		for(int n = 0; n < Nphases; n++)
		{
			Fractions(i,j,k)({n}) = 0.0;
		}
		for (auto it = Fields(i,j,k).cbegin();
				  it != Fields(i,j,k).cend(); ++it)
		{
			int pIndex = FieldsStatistics[it->index].Phase;
			Fractions(i, j, k)({pIndex}) += it->value;
		}
	OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::CalculateLaplacians(void)
{
	double dx_2 = 1.0 / (dx*dx);
	const int offset = Fields.Bcells() - 1;
	Storage3D<   Node, 0 >  FieldsCopy;
	FieldsCopy.Allocate(Fields.sizeX(), Fields.sizeY(), Fields.sizeZ(), Fields.Bcells());
	FieldsCopy = Fields;

	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, offset, )
	{
		if (Fields(i, j, k).flag)
		{
			for (auto it = Fields(i,j,k).begin(); it != Fields(i, j, k).end(); ++it)
			{
				it->laplacian = LaplacianStencil27[1][1][1] * (it->value)*dx_2;
			}
			for (int ii = -1; ii <= +1; ++ii)
			for (int jj = -1; jj <= +1; ++jj)
			for (int kk = -1; kk <= +1; ++kk)
			if (ii != 0 or jj != 0 or kk != 0)
			{
				for (auto it = FieldsCopy(i + ii, j + jj, k + kk).cbegin();
					it != FieldsCopy(i + ii, j + jj, k + kk).cend(); ++it)
					if (it->value != 0.0)
					{
						double value2 = LaplacianStencil27[ii + 1][jj + 1][kk + 1] * (it->value)*dx_2;
						Fields(i, j, k).add2(it->index, value2);
					}
			}
		}
	}
	OMP_PARALLEL_STORAGE_LOOP_END
}

NodeV PhaseField::Gradients(const int i, const int j, const int k) const
{
    double GWeights[3] = {-0.5/dx, 0.0, 0.5/dx};

    NodeV locGradients;
    for (int ii = -1; ii <= +1; ii += 2)
    {
        for (auto it = Fields(i+ii,j,k).cbegin();
                  it < Fields(i+ii,j,k).cend(); ++it)
        {
            locGradients.add_X(it->index, GWeights[ii+1]*(it->value));
        }
        for (auto it = Fields(i,j+ii,k).cbegin();
                  it < Fields(i,j+ii,k).cend(); ++it)
        {
            locGradients.add_Y(it->index, GWeights[ii+1]*(it->value));
        }
        for (auto it = Fields(i,j,k+ii).cbegin();
                  it < Fields(i,j,k+ii).cend(); ++it)
        {
            locGradients.add_Z(it->index, GWeights[ii+1]*(it->value));
        }
    }
    return locGradients;
}

NodeV PhaseField::Normals(const int i, const int j, const int k) const
{
    NodeV locNormal;
    NodeV locGradients = Gradients(i,j,k);

    for (auto alpha = Fields(i,j,k).cbegin();
              alpha != Fields(i,j,k).cend(); ++alpha)
    for (auto  beta = locGradients.cbegin();
               beta != locGradients.cend(); ++beta)
    if(alpha->index != beta->index)
    {

        locNormal.add(alpha->index, beta->index, alpha->value*beta->X,
                                                 alpha->value*beta->Y,
                                                 alpha->value*beta->Z);
    }

    for (auto alpha = locGradients.cbegin();
              alpha != locGradients.cend(); ++alpha)
    for (auto beta = Fields(i,j,k).cbegin();
              beta != Fields(i,j,k).cend(); ++beta)
    if(alpha->index != beta->index)
    {
        locNormal.add(alpha->index, beta->index, -alpha->X*beta->value,
                                                 -alpha->Y*beta->value,
                                                 -alpha->Z*beta->value);
    }

    for (auto alpha = locNormal.begin();
              alpha != locNormal.end(); ++alpha)
    {
        double Norm = sqrt(alpha->X*alpha->X +
                           alpha->Y*alpha->Y +
                           alpha->Z*alpha->Z);

        if (Norm > DBL_EPSILON)

        {
            Norm = 1.0/Norm;
            alpha->X *= Norm;
            alpha->Y *= Norm;
            alpha->Z *= Norm;
        }
        else
        {
            alpha->X = 0.0;
            alpha->Y = 0.0;
            alpha->Z = 0.0;
        }
    }
    return locNormal;
}

NodeV PhaseField::NormalsPhase(const int i, const int j, const int k) const
{
    NodeV locNormals = Gradients(i,j,k);

    for (auto alpha = locNormals.begin();
              alpha != locNormals.end(); ++alpha)
    {
        double Norm = sqrt(alpha->X*alpha->X +
                           alpha->Y*alpha->Y +
                           alpha->Z*alpha->Z);

        if (Norm > DBL_EPSILON)

        {
            Norm = 1.0/Norm;
            alpha->X *= Norm;
            alpha->Y *= Norm;
            alpha->Z *= Norm;
        }
        else
        {
            alpha->X = 0.0;
            alpha->Y = 0.0;
            alpha->Z = 0.0;
        }
    }
    return locNormals;
}

void PhaseField::CalculateVolumes(void)
{
    CalculateVolumesStepOne();
    CalculateVolumesStepTwo();
}

void PhaseField::CalculateVolumesStepOne(void)
{
    const size_t size = FieldsStatistics.size();
    for(size_t idx = 0; idx < size; idx++)
    {
        FieldsStatistics[idx].Volume = 0.0;
    }

    #ifdef _OPENMP
    #pragma omp parallel
    #endif
    {
        std::vector<double> Volume (size,0.0);
        OMP_PARALLEL_REDUCTION_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
        {
            if(Interface(i,j,k))
            {
                for(auto it = Fields(i,j,k).cbegin();
                         it < Fields(i,j,k).cend(); ++it)
                {
                    Volume[it->index] += it->value;
                }
            }
            else
            {
                const size_t idx = Fields(i,j,k).front().index; 
                Volume[idx] += 1.0;
            }
        }
        OMP_PARALLEL_REDUCTION_STORAGE_LOOP_END

        for(size_t idx = 0; idx < size; idx++)
        {
            //NOTE: atomic lead to arithmetic error in my code (Raphael Schiedung)
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            FieldsStatistics[idx].Volume += Volume[idx];
        }
    }
}

void PhaseField::CalculateVolumesStepTwo(void)
{
    for(unsigned int idx = 0; idx < FieldsStatistics.size(); idx++)
    {
        FieldsStatistics[idx].MAXVolume =
            max(FieldsStatistics[idx].Volume, FieldsStatistics[idx].MAXVolume);
    }
    int NumberOfNuclei = 0;

    for(unsigned int idx = 0; idx < FieldsStatistics.size(); idx++)
    if(FieldsStatistics[idx].Exist && FieldsStatistics[idx].Volume <= 0.0)
    {
        FieldsStatistics[idx].Exist  = 0;
        FieldsStatistics[idx].Stage  = 0;
        FieldsStatistics[idx].Volume = 0.0;
    }
    else if(FieldsStatistics[idx].Volume > 0.0)
    {
        FieldsStatistics[idx].Exist = 1;

        if((FieldsStatistics[idx].Stage == 1) &&
           (FieldsStatistics[idx].Volume > RefVolume))
        {
            FieldsStatistics[idx].Stage = 0;
        }
        NumberOfNuclei += FieldsStatistics[idx].Stage;
    }
    NucleationPresent = (NumberOfNuclei > 0);
}

void PhaseField::SetFlags(void)
{
    const int offset = Fields.Bcells()-1;

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, offset,)
    {
        if(Fields(i,j,k).flag == 2)
        {
            for(int ii = -1; ii <= +1; ++ii)
            for(int jj = -1; jj <= +1; ++jj)
            for(int kk = -1; kk <= +1; ++kk)
            if(ii != 0 or jj != 0 or kk != 0)
            if(!(Fields(i+ii, j+jj, k+kk).flag))
            {
                Fields(i+ii, j+jj, k+kk).flag = 1;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::Finalize(const BoundaryConditions& BC, bool finalize)
{
    if(finalize)
    {
        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
        {
            if(Fields(i,j,k).flag)
            {
                Fields(i,j,k).finalize();
            }
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    SetBoundaryConditions(BC);
    SetFlags();
    CalculateLaplacians();
    CalculateVolumes();
	CalculateFractions();
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
/*
void PhaseField::NormalizeIncrements(const BoundaryConditions& BC, const double dt)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    if (Fields(i,j,k).flag)
    {
        if(FieldsDot(i,j,k).size() == 1)
        {
            double limit = 1.0;
            int    indexA     = FieldsDot(i,j,k).front().indexA;
            double increment  = FieldsDot(i,j,k).front().value * dt;
            double PFvalue    = Fields(i,j,k)[indexA];
            double newPFvalue = PFvalue + increment;
            if(newPFvalue < 0.0 and PFvalue > 0.0)
            {
                limit = min(limit, -newPFvalue/increment);
                FieldsDot(i,j,k).front().value *= limit;
            }
            else if(newPFvalue > 1.0 and PFvalue < 1.0)
            {
                limit = min(limit,(1.0 - newPFvalue)/increment);
                FieldsDot(i,j,k).front().value *= limit;
            }
            else
            {
                FieldsDot(i,j,k).clear();
            }
        }
        else
        {
            for(auto alpha = Fields(i,j,k).cbegin();
                     alpha != Fields(i,j,k).cend(); ++alpha)
            {
                // Set clearly corrupted sets to zero

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
                // Remove zero-sets from the storage

                if(alpha->value == 0.0)
                {
                    alpha = FieldsDot(i,j,k).erase(alpha);
                }
                else
                {
                    ++alpha;
                }
            }

            // Limit increments! This is done in a while loop, to acknowledge all
            // existing pair-contributions.

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
                    // Collect increments

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
                    // Calculate limits

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
                    // Limit increments

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
                // cleaning up values which are too small

                if(fabs(alpha->value) < DBL_EPSILON)
                {
                    alpha = FieldsDot(i,j,k).erase(alpha);
                }
                else
                {
                    ++alpha;
                }
            }

            // Skip limiting loop, if no convergence after 24 iterations

            if (number_of_iterations > 24)
            LimitingNeeded = false;

            // End plausibility check

            Node locField = Fields(i,j,k);
            for(auto psi = FieldsDot(i,j,k).cbegin();
                     psi < FieldsDot(i,j,k).cend(); ++psi)
            if(psi->value != 0.0)
            {
                locField.add(psi->indexA,  psi->value * dt);
                locField.add(psi->indexB, -psi->value * dt);
            }
            for(auto it = locField.cbegin();
                     it < locField.cend(); ++it)
            if(it->value < -DBL_EPSILON or it->value > 1.0+DBL_EPSILON)
            {
                double oldFields = Fields(i,j,k).get(it->index);
                string msg = "Normalizing of phase field increments failed in point ("
                           + to_string(i) + "," + to_string(j) + "," + to_string(k)
                           + "). " + to_string(Fields(i,j,k).size())
                           + " fields present. Grain "
                           + to_string(it->index) + " with a fields-value of "
                           + to_string(oldFields)
                           + " is incremented by "
                           + to_string(it->value-oldFields)
                           + ", which results in a fields-value of "
                           + to_string(it->value)
                           + ". This will result in undefined behavior!";
                Info::WriteWarning(msg,this->thisclassname,"NormalizeIncrements");
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    SetIncrementsBoundaryConditions(BC);
}*/

void PhaseField::SetBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(Fields);
    BC.SetY(Fields);
    BC.SetZ(Fields);
}
void PhaseField::SetIncrementsBoundaryConditions(const BoundaryConditions& BC)
{
    BC.SetX(FieldsDot);
    BC.SetY(FieldsDot);
    BC.SetZ(FieldsDot);
}
void PhaseField::PrintPointStatistics(const int x, const int y, const int z) const
{
    cout << "Point: " << x << " " << y << " " << z << endl;
    cout << "Phase Fields: " << endl;
    cout << "Index: ";
    for (auto beta = Fields(x,y,z).cbegin();
              beta != Fields(x,y,z).cend();  ++beta)
    {
        cout << setw(12) << setfill(' ') << beta->index << " ";
    }
    cout << endl;
    cout << "Value: ";
    for (auto beta = Fields(x,y,z).cbegin();
              beta != Fields(x,y,z).cend();  ++beta)
    {
        cout.precision(6);
        cout << setw(12) << setfill(' ') << beta->value << " ";
    }
    cout << endl;
    cout << "Phase Fractions: " << endl;
    cout << "Index: ";
    for (auto alpha = 0; alpha != Nphases; alpha++)
    {
        cout << setw(12) << setfill(' ') << alpha << " ";
    }
    cout << endl;
    cout << "Value: ";

    auto locFractions = Fractions(x,y,z);
    for (auto alpha = 0; alpha != Nphases; alpha++)
    {
        cout.precision(6);
        cout << setw(12) << setfill(' ') << locFractions({alpha}) << " ";
    }
    cout << endl;
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

void PhaseField::WriteVTKData(std::stringstream& buffer,
        const bool CurvatureOutput) const
{
    // Interface
    buffer << "<DataArray type = \"Float64\" Name = \"Interfaces\""
           << " NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        double sum = 0.5;

        if(Interface(i,j,k))
        for (auto alpha = Fields(i,j,k).cbegin();
                  alpha != Fields(i,j,k).cend() - 1; ++alpha)
        for (auto  beta = alpha + 1;
                   beta < Fields(i,j,k).cend(); ++beta)
        {
            sum -= alpha->value*beta->value;
        }
        buffer << 1.0/(sum*2.0) << "\n";
    }
    buffer << "</DataArray>" << endl;
    // Flags
    buffer << "<DataArray type = \"Float64\" Name = \"Flags\""
           << " NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        if(Fields(i,j,k).flag)
        {
            buffer << Fields(i,j,k).flag << "\n";
        }
        else
        {
            buffer << 0 << "\n";
        }
    }
    buffer << "</DataArray>" << endl;
    // Phase-Fields
    buffer << "<DataArray type = \"Float64\" Name = \"PhaseFields\""
           << " NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        double sum = 0.0;
        for (auto beta = Fields(i,j,k).cbegin();
                  beta < Fields(i,j,k).cend();  ++beta)
        {
            sum += beta->value*beta->index;
        }
        buffer << sum << "\n";
    }
    buffer << "</DataArray>" << endl;
    // Phase Fractions
    for(int n = 0; n < Nphases; n++)
    {
        buffer << "<DataArray type = \"Float64\" Name = \"PhaseFraction_" << n << "\""
               << " NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for (int k = 0; k < Nz; k++)
        for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++)
        {
            buffer << Fractions(i,j,k)({n}) << "\n";
        }
        buffer << "</DataArray>" << endl;
    }
    // Local number of phase fields
    buffer << "<DataArray type = \"Float64\" Name = \"Junctions\""
           << " NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for (int k = 0; k < Nz; k++)
    for (int j = 0; j < Ny; j++)
    for (int i = 0; i < Nx; i++)
    {
        double sum = 0.0;

        if(Interface(i,j,k))
        {
            for (auto alpha = Fields(i,j,k).cbegin();
                      alpha != Fields(i,j,k).cend();  ++alpha)
            if(alpha->value != 0.0)
            {
                sum += 1.0;
            }
        }
        else
        {
            sum = 1.0;
        }
        buffer << sum << "\n";
    }
    buffer << "</DataArray>" << endl;
    buffer << "<DataArray type = \"Float64\" Name = \"Variants\""
           << " NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        int locVariant = 0;
        double locValue = 0.0;
        for (auto beta = Fields(i,j,k).cbegin();
                  beta != Fields(i,j,k).cend(); ++beta)
        {
            if(beta->value > locValue)
            {
                locValue = beta->value;
                locVariant = FieldsStatistics[beta->index].Variant;
            }
        }
        buffer << locVariant << "\n";
    }
    buffer << "</DataArray>" << endl;
    if(CurvatureOutput)
    {
        // Write extended output.
        // If the extended output is used the curvature of each phase will
        // be written in to the vtk file.
        for (int phase = 0; phase < Nphases; phase++)
        {
            buffer << "<DataArray type = \"Float64\" Name = "
                <<"\"Mean curvature(" + to_string(phase) + ")\""
                << " NumberOfComponents=\"1\" format=\"ascii\">" << endl;

            for (int k = 0; k < Nz; k++)
            for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
            if (Interface(i,j,k))
            {
                Node   locKappa1; // principle curvature of each phase field
                Node   locKappa2; // principle curvature of each phase field
                double phaseKappa1 = 0; // principle curvature phase
                double phaseKappa2 = 0; // principle curvature phase

                // Calculate curvature of all phase fields
                tie(locKappa1, locKappa2) = PrincipleCurvatures(i,j,k);

                // Calculate curvature of desired phase
                for (auto it = locKappa1.cbegin();
                        it < locKappa1.cend(); it++)
                {
                    int pIndex = FieldsStatistics[it->index].Phase;
                    if (pIndex == phase) phaseKappa1 += it->value;
                }
                for (auto it = locKappa2.cbegin();
                        it < locKappa2.cend(); it++)
                {
                    int pIndex = FieldsStatistics[it->index].Phase;
                    if (pIndex == phase) phaseKappa2 += it->value;
                }
                buffer << 0.5 * (phaseKappa1 + phaseKappa2) << "\n";
            }
            else buffer << 0 << "\n";
            buffer << "</DataArray>" << endl;
            buffer << "<DataArray type = \"Float64\" Name = "
                   <<"\"Principle curvature(" + to_string(phase) + ")\""
                   << " NumberOfComponents=\"2\" format=\"ascii\">" << endl;

            for (int k = 0; k < Nz; k++)
            for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
            if (Interface(i,j,k))
            {
                Node   locKappa1; // principle curvature of each phase field
                Node   locKappa2; // principle curvature of each phase field
                double phaseKappa1 = 0; // principle curvature phase
                double phaseKappa2 = 0; // principle curvature phase

                tie(locKappa1, locKappa2) = PrincipleCurvatures(i,j,k);
                for (auto it = locKappa1.cbegin();
                          it < locKappa1.cend(); it++)
                {
                    int pIndex = FieldsStatistics[it->index].Phase;
                    if (pIndex == phase) phaseKappa1 += it->value;
                }
                for (auto it = locKappa2.cbegin();
                          it < locKappa2.cend(); it++)
                {
                    int pIndex = FieldsStatistics[it->index].Phase;
                    if (pIndex == phase) phaseKappa2 += it->value;
                }
                buffer << phaseKappa1 << " " << phaseKappa2 << "\n";
            }
            else buffer << 0 << " " << 0 << "\n";
            buffer << "</DataArray>" << endl;
        }
    }
}


void PhaseField::WriteDistortedVTK(int tStep, ElasticitySolverSpectral& ES,
        double scale)
{
    stringstream outbufer;

    outbufer << "# vtk DataFile Version 3.0\n";
    outbufer << "PhaseField with distortions\n";
    outbufer << "ASCII\n";
    outbufer << "DATASET STRUCTURED_GRID\n";
    outbufer << "DIMENSIONS " << Nx << " " << Ny << " " << Nz << "\n";
    outbufer << "POINTS " <<  Nx*Ny*Nz << " double\n";

    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        outbufer << i + scale * ES.rlU[0][k + Nz*(j + Ny*i)]/dx << " "
                 << j + scale * ES.rlU[1][k + Nz*(j + Ny*i)]/dx << " "
                 << k + scale * ES.rlU[2][k + Nz*(j + Ny*i)]/dx << "\n";
    }
    outbufer << "\n";
    outbufer << "POINT_DATA " << Nx*Ny*Nz << "\n";

    outbufer << "SCALARS Interfaces double 1" << "\n";
    outbufer << "LOOKUP_TABLE default" << "\n";
    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        double sum = 0.5;
        if(Interface(i,j,k))
        for (auto alpha = Fields(i,j,k).cbegin();
                  alpha != Fields(i,j,k).cend() - 1; ++alpha)
        for (auto  beta = alpha + 1;
                   beta < Fields(i,j,k).cend(); ++beta)
        {
            sum -= alpha->value*beta->value;
        }
        outbufer << 1.0/(sum*2.0) << "\n";
    }
    outbufer << "SCALARS PhaseFields double 1" << "\n";
    outbufer << "LOOKUP_TABLE default" << "\n";
    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        double sum = 0.0;
        for (auto beta = Fields(i,j,k).cbegin();
                  beta < Fields(i,j,k).cend();  ++beta)
        {
            sum += beta->value*beta->index;
        }
        outbufer << sum << "\n";
    }
    for (int n = 0; n < Nphases; n++)
    {
        outbufer << "SCALARS PhaseFraction_" << n << " double 1" << "\n";
        outbufer << "LOOKUP_TABLE default" << "\n";
        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            outbufer << Fractions(i,j,k)({n}) << "\n";
        }
    }
    outbufer << "SCALARS ColorScale double 1" << "\n";
    outbufer << "LOOKUP_TABLE default" << "\n";
    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        int locPF = 0;
        double locVal = 0.0;
        for (auto beta = Fields(i,j,k).cbegin();
                  beta < Fields(i,j,k).cend();  ++beta)
        {
            if(beta->value > locVal)
            {
                locVal = beta->value;
                locPF = beta->index;
            }
        }
        outbufer << locPF << "\n";
    }
    outbufer << "SCALARS Variants double 1" << "\n";
    outbufer << "LOOKUP_TABLE default" << "\n";
    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        int locVariant = 0;
        double locValue = 0.0;
        for (auto beta = Fields(i,j,k).cbegin();
                  beta != Fields(i,j,k).cend(); ++beta)
        {
            if(beta->value > locValue)
            {
                locValue = beta->value;
                locVariant = FieldsStatistics[beta->index].Variant;
            }
        }
        outbufer << locVariant << "\n";
    }
    outbufer << "SCALARS Flags double 1" << "\n";
    outbufer << "LOOKUP_TABLE default" << "\n";
    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        outbufer << Fields(i,j,k).flag << "\n";
    }

    string FileName = UserInterface::MakeFileName(VTKDir,"DistPhaseField_",
            tStep, ".vtk");

    ofstream vtk_file(FileName.c_str());
    vtk_file << outbufer.rdbuf();
    vtk_file.close();
}

void PhaseField::WriteLaplaciansVTK(int tStep, int PhiIndex)
{
    stringstream buffer;

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    std::vector<int> DataTypes {PDScalars};
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        buffer << "<DataArray type = \"Float64\" "
               << "Name = \"Laplacian\" NumberOfComponents=\"1\" "
               << "format=\"ascii\">" << endl;
        for (int k = 0; k < Nz; k++)
        for (int j = 0; j < Ny; j++)
        for (int i = 0; i < Nx; i++)
        {
            if(Fields(i,j,k).flag)
            {
                buffer << Fields(i,j,k).get2(PhiIndex) << endl;
            }
            else
            {
                buffer << 0.0 << endl;
            }
        }
        buffer << "</DataArray>" << endl;
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "Laplacians_" + to_string(PhiIndex), tStep);
}

void PhaseField::WriteIndividualPhaseFieldValuesVTK(const int tStep,
        const std::initializer_list<int> FieldIndices) const
{
    stringstream buffer;

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    std::vector<int> DataTypes {PDScalars};
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        for (auto IteratorIndex = FieldIndices.begin();
                  IteratorIndex != FieldIndices.end(); IteratorIndex++)
        if(FieldsStatistics[*IteratorIndex].Exist == true)
        {
            buffer << "<DataArray type = \"Float64\" Name = \"FieldValue ("
                << *IteratorIndex << ")\" NumberOfComponents=\"1\" "
                << "format=\"ascii\">" << endl;
            for (int k = 0; k < Nz; k++)
            for (int j = 0; j < Ny; j++)
            for (int i = 0; i < Nx; i++)
            {
                double val = 0.0;

                if(Interface(i,j,k))
                {
                    for (auto alpha = Fields(i,j,k).cbegin();
                              alpha != Fields(i,j,k).cend(); ++alpha)
                    {
                        if(alpha->index == *IteratorIndex){val = alpha->value;};
                    }
                }
                else
                {
                    if(Fields(i,j,k).front().index == *IteratorIndex){val = 1.0;};
                }
                buffer << val << "\n";
            }
            buffer << "</DataArray>" << endl;
        }
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "PhaseFieldValues", tStep);
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

void PhaseField::Read(const BoundaryConditions& BC, int tStep)
{
    string FileName =
        UserInterface::MakeFileName(RawDataDir,"PF_", tStep, ".dat");

    Read(FileName, BC);
    FieldsStatistics.Read(tStep);
    Finalize(BC);
}

void PhaseField::Read(string FileName, const BoundaryConditions& BC)
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

    Finalize(BC);
    Info::WriteStandard(thisclassname, "Binary input loaded");
}

void PhaseField::WriteAverageVolume(const int tStep, const int PhaseIndex) const
{
    stringstream converter;
    converter << PhaseIndex;

    string FileName = string("AverageVolumeOfPhase_")
        + converter.str() + string(".txt");

    if(!tStep)
    {
        fstream tout(FileName.c_str(), ios::out);
        tout << "Time\tAvgVolume" << endl;
        tout.close();
    }

    double avgVol = 0;
    double count = 0;
    for(unsigned int n = 0; n < FieldsStatistics.size(); n++)
    if(FieldsStatistics[n].Exist and FieldsStatistics[n].Phase == PhaseIndex)
    {
        avgVol += FieldsStatistics[n].Volume;
        count += 1.0;
    }
    fstream out(FileName.c_str(), ios::out | ios::app);
    if(count)
    {
        out << tStep << "\t" << avgVol/count << endl;
    }
    else
    {
        out << tStep << "\t" << 0.0 << endl;
    }
    out.close();
}

void PhaseField::MoveFrame(const int dx, const int dy, const int dz,
                           const BoundaryConditions& BC)
{
    int xBeg = (dx >= 0) + (dx < 0)*(Nx) - 1;
    int xEnd = (dx >= 0)*(Nx) + (dx < 0) - 1;
    int xInc = 1 - 2*(dx < 0);

    int yBeg = (dy >= 0) + (dy < 0)*(Ny) - 1;
    int yEnd = (dy >= 0)*(Ny) + (dy < 0) - 1;
    int yInc = 1 - 2*(dy < 0);

    int zBeg = (dz >= 0) + (dz < 0)*(Nz) - 1;
    int zEnd = (dz >= 0)*(Nz) + (dz < 0) - 1;
    int zInc = 1 - 2*(dz < 0);

    for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
    {
        Fields(i, j, k) = Fields(i+dx,j+dy,k+dz);
    }
    Finalize(BC);
    Info::WriteStandard(thisclassname, "Frame moved");
}


void PhaseField::ConsumePlane(const int dx, const int dy, const int dz,
                           const int x, const int y, const int z,
                           const BoundaryConditions& BC)
{
    int xBeg = (dx >= 0) + (dx < 0)*(Nx) - 1;
    int xEnd = (dx >= 0)*(Nx) + (dx < 0) - 1;
    int xInc = 1 - 2*(dx < 0);

    int yBeg = (dy >= 0) + (dy < 0)*(Ny) - 1;
    int yEnd = (dy >= 0)*(Ny) + (dy < 0) - 1;
    int yInc = 1 - 2*(dy < 0);

    int zBeg = (dz >= 0) + (dz < 0)*(Nz) - 1;
    int zEnd = (dz >= 0)*(Nz) + (dz < 0) - 1;
    int zInc = 1 - 2*(dz < 0);

    for(int i = xBeg; ((dx >= 0) and (i <= xEnd)) or ((dx < 0) and (i >= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j <= yEnd)) or ((dy < 0) and (j >= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k <= zEnd)) or ((dz < 0) and (k >= zEnd)); k += zInc)
    if((i-x)*dx+(j-y)*dy+(k-z)*dz >= 0)
    {
        Fields(i, j, k) = Fields(i+dx,j+dy,k+dz);
    }
    xBeg = (dx >= 0)*(Nx) + (dx < 0) - 1;
    xEnd = (dx >= 0) + (dx < 0)*(Nx) - 1;
    xInc = 2*(dx < 0) - 1;

    yBeg = (dy >= 0)*(Ny) + (dy < 0) - 1;
    yEnd = (dy >= 0) + (dy < 0)*(Ny) - 1;
    yInc = 2*(dy < 0) - 1;

    zBeg = (dz >= 0)*(Nz) + (dz < 0) - 1;
    zEnd = (dz >= 0) + (dz < 0)*(Nz) - 1;
    zInc = 2*(dz < 0) - 1;

    for(int i = xBeg; ((dx >= 0) and (i >= xEnd)) or ((dx < 0) and (i <= xEnd)); i += xInc)
    for(int j = yBeg; ((dy >= 0) and (j >= yEnd)) or ((dy < 0) and (j <= yEnd)); j += yInc)
    for(int k = zBeg; ((dz >= 0) and (k >= zEnd)) or ((dz < 0) and (k <= zEnd)); k += zInc)
    if((i-x)*dx+(j-y)*dy+(k-z)*dz < 0)
    {
        Fields(i, j, k) = Fields(i-dx,j-dy,k-dz);
    }
    Finalize(BC);
    Info::WriteStandard(thisclassname, "Plane consumed");
}

//void PhaseField::CalculateGrainboundaryPotentialContrib(ThermodynamicFunctions& TF, Settings& locSettings)
//{
//    int Ncomp = locSettings.Ncomp;
//    STORAGE_LOOP_BEGIN(i,j,k,Fields,0)
//    if(Interface(i,j,k))
//    for(auto alpha = Fields(i, j, k).cbegin();
//                     alpha != Fields(i, j, k).cend(); ++alpha)
//    for(int comp = 0; comp < Ncomp; comp++)
//    {
//        double contributionfactor = 5E-1;
//        int index = alpha->index;
//        int pIndex = FieldsStatistics[index].Phase;
//        TF.ChemicalPotential(i,j,k)({pIndex,comp}) -= fabs(
//                                     (1.0-Fractions(i,j,k)({pIndex}))
//                                     *TF.ChemicalPotential(i,j,k)({pIndex,comp})
//                                     *contributionfactor);
//    }
//    STORAGE_LOOP_END
//}

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

void PhaseField::PrintVolumeFractions(std::vector<std::string> Names)
{
    vector<double> tPhaseFrac(Nphases);
    double totalVolume = double(Nx * Ny * Nz);
    for(unsigned int idx = 0; idx < FieldsStatistics.size(); idx++)
    {
        tPhaseFrac[FieldsStatistics[idx].Phase] += FieldsStatistics[idx].Volume;
    }
    Info::WriteSimple("Phase fractions:");
    for(unsigned int idx = 0; idx < tPhaseFrac.size(); idx++)
    {
        Info::WriteStandardNarrow("Phase " + std::to_string(idx) + " (" + Names[idx] + ")",
                std::to_string( (100.0*tPhaseFrac[idx])/totalVolume ) + " %" );
    }
}

void PhaseField::PrintVolumeFractions(ChemicalProperties& CP)
{
    vector<double> tPhaseFrac(Nphases);
    double totalVolume = double(Nx * Ny * Nz);
    for(unsigned int idx = 0; idx < FieldsStatistics.size(); idx++)
    {
        tPhaseFrac[FieldsStatistics[idx].Phase] += FieldsStatistics[idx].Volume;
    }
    Info::WriteSimple("Phase fractions:");
    for(unsigned int idx = 0; idx < tPhaseFrac.size(); idx++)
    {
        Info::WriteStandardNarrow("Phase " + std::to_string(idx) + " (" + CP.Phase[idx].Name +")",
                std::to_string( (100.0*tPhaseFrac[idx])/totalVolume ) + " %" );
    }
}

void PhaseField::Remesh(int newNx, int newNy, int newNz,
        const BoundaryConditions& BC)
{
#ifdef _OPENMP
    destroy_omp_locks();
#endif
    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    Fields.Remesh(Nx, Ny, Nz);
    FieldsDot.Reallocate(Nx, Ny, Nz);
	Fractions.Reallocate(Nx, Ny, Nz);
#ifdef _OPENMP
    omp_locks.Reallocate(Nx,Ny,Nz);
    init_omp_locks();
#endif

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        Fields(i,j,k).flag = 2*(Fields(i,j,k).size() > 1);
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Finalize(BC);

    RefVolume = Pi;

    if(Nx < iWidth) RefVolume *= (Nx+1)/2.0;
    else RefVolume *= 1.1*iWidth;
    if(Ny < iWidth) RefVolume *= (Ny+1)/2.0;
    else RefVolume *= 1.1*iWidth;
    if(Nz < iWidth) RefVolume *= (Nz+1)/2.0;
    else RefVolume *= 1.1*iWidth;

    Info::WriteStandard(thisclassname, "Remeshed");
}

int PhaseField::PlantGrainNucleus(int PhaseIndex, int x, int y, int z)
{
    int locIndex = FieldsStatistics.add_nucleus(PhaseIndex);
    FieldsStatistics[locIndex].Rcm[0] = x;
    FieldsStatistics[locIndex].Rcm[1] = y;
    FieldsStatistics[locIndex].Rcm[2] = z;
    
    NucleationPresent = true;

    int Range = 1;
    for(int i = -Range; i <= Range; i++)
    for(int j = -Range; j <= Range; j++)
    for(int k = -Range; k <= Range; k++)
    if(i*i + j*j + k*k <= 2*Range*Range)
    {
        Fields(x+i, y+j, z+k).flag = max(1, Fields(x+i, y+j, z+k).flag);
    }
    Fields(x, y, z).set(locIndex, 0.0);
    return locIndex;
}

int PhaseField::AddGrainInfo(int PhaseIndex)
{
    int locIndex = FieldsStatistics.add_grain(PhaseIndex);

    return locIndex;
}

PhaseField& PhaseField::operator= (const PhaseField& rhs)
{
    // protect against invalid self-assignment and copy of unitialized object
    if (this != &rhs and rhs.thisclassname == "PhaseField")
    {
        thisclassname = rhs.thisclassname;
        NucleationPresent = rhs.NucleationPresent;
        Nx = rhs.Nx;
        Ny = rhs.Ny;
        Nz = rhs.Nz;
        dx  = rhs.dx;
        Eta = rhs.Eta;
        Nphases = rhs.Nphases;
        iWidth = rhs.iWidth;
        RefVolume = rhs.RefVolume;

        if (FieldsStatistics.size() != rhs.FieldsStatistics.size())
        {
            FieldsStatistics.Allocate(rhs.FieldsStatistics.size());
        }

        if (Fields.IsNotAllocated())
        {
            FieldsStatistics.Allocate(rhs.FieldsStatistics.size());

            Fields.Allocate(rhs.Nx, rhs.Ny, rhs.Nz, rhs.Fields.Bcells());
            FieldsDot.Allocate(rhs.Nx, rhs.Ny, rhs.Nz, {rhs.Nphases},
                    rhs.FieldsDot.Bcells());
        }
        else if (not Fields.IsSize(rhs.Nx, rhs.Ny, rhs.Nz))
        {
            Fields.Reallocate(rhs.Nx, rhs.Ny, rhs.Nz);
            FieldsDot.Reallocate(rhs.Nx, rhs.Ny, rhs.Nz);
        }

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
        {
            Fields(i,j,k) = rhs.Fields(i,j,k);
            FieldsDot(i,j,k) = rhs.FieldsDot(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END

        FieldsStatistics = rhs.FieldsStatistics;

#ifdef _OPENMP
        omp_locks.Allocate(Nx, Ny, Nz, Fields.Bcells());
        init_omp_locks();
#endif
    }
    return *this;
}

std::vector<int> PhaseField::GetPresentPhaseFields() const
{
    std::vector<int> indeces;
    for(unsigned int n = 0; n < FieldsStatistics.size(); n++)
    if(FieldsStatistics[n].Exist)
    {
        indeces.push_back(n);
    }

    return indeces;
}

std::vector<int> PhaseField::ReturnVicinityPhaseFields(const int i,
                                           const int j, const int k) const
{
    // Analyze vicinity
    std::vector<int> tempPFindex;
    tempPFindex.push_back (Fields(i,j,k).front().index);
    for(int ii = -1; ii <= 1; ii++)
    for(int jj = -1; jj <= 1; jj++)
    for(int kk = -1; kk <= 1; kk++)
    {
        for(auto alpha = Fields(i+ii,j+jj,k+kk).cbegin();
                  alpha != Fields(i+ii,j+jj,k+kk).cend(); ++alpha)
        {
            if(alpha -> index != Fields(i,j,k).front().index)
            {
                tempPFindex.push_back (alpha->index);
            }
        }
    }

    // Erase double entries
    std::sort(tempPFindex.begin(), tempPFindex.end() );
    tempPFindex.erase( std::unique( tempPFindex.begin(),
                tempPFindex.end() ), tempPFindex.end() );

    return tempPFindex;
}

void PhaseField::CombinePhaseFields(const int PhaseIndex)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
    {
        double phasevalue = 0.0;
        bool ChangedPF = false;
        if (Interface(i,j,k))
        {
            for (auto it = Fields(i,j,k).begin();
                      it < Fields(i,j,k).end(); ++it)
            {
                if(FieldsStatistics[it->index].Phase and it->index != PhaseIndex)
                {
                    it -> index = PhaseIndex;
                    phasevalue += it -> value;
                    it -> value = 0.0;
                    ChangedPF = true;
                }
            }
        }
        else
        {
            if(FieldsStatistics[Fields(i,j,k).front().index].Phase and
                    Fields(i,j,k).front().index != PhaseIndex)
            {
                Fields(i,j,k).begin()->index = PhaseIndex;
                phasevalue += Fields(i,j,k).begin()->value;
                Fields(i,j,k).begin()->value = 0.0;
                Fields(i,j,k).front().index = PhaseIndex;
                ChangedPF = true;
            }
        }

        if (ChangedPF == true) Fields(i,j,k).set(PhaseIndex,
                                Fields(i,j,k).get(PhaseIndex) + phasevalue);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void PhaseField::SelectiveCombinePhaseFields(BoundaryConditions& BC, const int SourcePhaseIndex,
        const int TargetPhaseIndex)
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,Fields.Bcells(),)
    {
        double phasevalue = 0.0;
        bool ChangedPF = false;
        if (Interface(i,j,k))
        {
            for (auto it = Fields(i,j,k).begin();
                      it < Fields(i,j,k).end(); ++it)
            {
                if(it->index == SourcePhaseIndex)
                {
                    it -> index = TargetPhaseIndex;
                    phasevalue += it -> value;
                    it -> value = 0.0;
                    ChangedPF = true;
                }
            }
        }
        else
        {
            if(Fields(i,j,k).front().index == SourcePhaseIndex)
            {
                Fields(i,j,k).begin()->index = TargetPhaseIndex;
                phasevalue += Fields(i,j,k).begin()->value;
                Fields(i,j,k).begin()->value = 0.0;
                Fields(i,j,k).front().index = TargetPhaseIndex;
                ChangedPF = true;
                Fields(i,j,k).flag = 1;
            }
        }

        if (ChangedPF == true) Fields(i,j,k).set(TargetPhaseIndex,
                                Fields(i,j,k).get(TargetPhaseIndex) + phasevalue);
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    Finalize(BC);
}

pair<Node, Node> PhaseField::PrincipleCurvatures(const int i, const int j,
        const int k) const
{
    Node kappa1;
    Node kappa2;

    kappa1.clear();
    kappa2.clear();

    const double GWeights[3] = {-0.5/dx, 0.0, 0.5/dx};

    // Check if neighbour cell are also in the interface
    if (Interface(i,j,k))
    {
        // Calculate gradients of phase normal fields (Jacobian matrix)
        Tensor<Node, 2 > NormalGradient;
        NormalGradient.Allocate({3,3});
        NormalGradient.set_to_zero();

        // Calculate gradients of normals
        for (int ii = -1; ii <= +1; ii += 2)
        {
            const double GWeight = GWeights[ii+1];

            NodeV locNormalsX = NormalsPhase(i+ii,j,k);
            NodeV locNormalsY = NormalsPhase(i,j+ii,k);
            NodeV locNormalsZ = NormalsPhase(i,j,k+ii);

            for (auto it = locNormalsX.cbegin(); it < locNormalsX.cend(); ++it)
            {
                NormalGradient({0,0}).add(it->index, GWeight*(it->X));
                NormalGradient({1,0}).add(it->index, GWeight*(it->Y));
                NormalGradient({2,0}).add(it->index, GWeight*(it->Z));
            }
            for (auto it = locNormalsY.cbegin(); it < locNormalsY.cend(); ++it)
            {
                NormalGradient({0,1}).add(it->index, GWeight*(it->X));
                NormalGradient({1,1}).add(it->index, GWeight*(it->Y));
                NormalGradient({2,1}).add(it->index, GWeight*(it->Z));
            }
            for (auto it = locNormalsZ.cbegin(); it < locNormalsZ.cend(); ++it)
            {
                NormalGradient({0,2}).add(it->index, GWeight*(it->X));
                NormalGradient({1,2}).add(it->index, GWeight*(it->Y));
                NormalGradient({2,2}).add(it->index, GWeight*(it->Z));
            }
        }

        // Calculate basis of tangent space at (i,j,k)
        // the basis of a spherical coordinate system at (i,j,k) will be
        // used here. The phase normal vector is in this case equal to the
        // radial normal vector of the spherical coordinate and the
        // remaining basis vectors will be the basis of the tangent space

        NodeV locNormals = NormalsPhase(i,j,k);
        for (auto it = locNormals.cbegin(); it < locNormals.cend(); ++it)
        {
            double phi   = atan2(it->Y,it->X);
            double theta = acos(it->Z);

            // Calculate basis vectors of tangent space
            dVector3 e_n;
            e_n.set_to_zero();
            e_n[0] = it->X;
            e_n[1] = it->Y;
            e_n[2] = it->Z;

            dVector3 e_theta;
            e_theta.set_to_zero();
            e_theta[0] = cos(theta) * cos(phi);
            e_theta[1] = cos(theta) * sin(phi);
            e_theta[2] =            - sin(theta);

            dVector3 e_phi;
            e_phi.set_to_zero();
            e_phi[0] = - sin(phi);
            e_phi[1] =   cos(phi);

            // Calculate projection and inclusion matrices
            double Projection [2][3];
            double Inclusion  [3][2];
            for (int m = 0; m < 3; ++m)
            {
                Projection[1][m] = e_phi  [m];
                Projection[0][m] = e_theta[m];

                Inclusion[m][1] = Projection[1][m];
                Inclusion[m][0] = Projection[0][m];
            }

            // Calculate local Weingarten map W
            complex<double> W[2][2];
            for (int l = 0; l < 2; ++l)
            for (int m = 0; m < 2; ++m)
            {
                W[l][m] = 0.0;
                for (int n = 0; n < 3; ++n)
                for (int o = 0; o < 3; ++o)
                {
                    W[l][m] -= Projection[l][n] *
                        NormalGradient({n,o}).get(it->index)
                        * Inclusion[o][m];
                }
            }

            // Calculate eigenvalues of local Weingarten map
            const double part1 = real((W[0][0] + W[1][1])/2.);
            const double W2    = real((W[0][0]-W[1][1]) * (W[0][0]-W[1][1]));
            const double part2 = real(0.5 * sqrt(W2) + 4.*W[0][1] * W[1][0]);

            const double locKappa1 = real(part1 - part2);
            const double locKappa2 = real(part1 + part2);

            // Add eigenvalues to Node storage
            kappa1.add(it->index, locKappa1);
            kappa2.add(it->index, locKappa2);
        }
    }

    return make_pair(kappa1, kappa2);
}

bool PhaseField::PhaseFieldPresent(const int i, const int j, const int k,
                                                        const int Index) const
{
    if(!Interface(i,j,k))
    {
        if (Fields(i,j,k).front().index == Index) return true;
    }
    else
    {
        for (auto alpha = Fields(i,j,k).cbegin();
                  alpha != Fields(i,j,k).cend(); ++alpha)
        {
            if (alpha->index == Index) return true;
        }
    }
    return false;
}

std::vector<int> PhaseField::GetMaxPhaseFieldOverlap(const int thPhase1,
        const int thPhase2)
{
    std::vector<int> maxOverlap;
    Matrix<int> Overlap;
    int numberOfGrains = FieldsStatistics.GrainStorage.size();
    Overlap.Allocate(numberOfGrains, numberOfGrains);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        if (Interface(i,j,k))
        {
            for (auto it = Fields(i,j,k).cbegin();
                      it != Fields(i,j,k).cend(); ++it)
            {
                int locIndex1 = it->index;
                int thPhaseIndex1 = FieldsStatistics[locIndex1].Phase;
                if (thPhaseIndex1 == thPhase1)
                {
                    for (auto jt = it+1;
                              jt != Fields(i,j,k).cend(); ++jt)
                    {
                        int locIndex2 = jt->index;
                        int thPhaseIndex2 = FieldsStatistics[locIndex2].Phase;
                        if (thPhaseIndex2 == thPhase2)
                        {
                            Overlap.add(locIndex1, locIndex2, 1);
                        }
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    int maxNx = -1;
    int maxNy = -1;
    for (int it = 0; it < numberOfGrains-1; it++)
    for (int jt = it+1; jt < numberOfGrains; jt++)
    {
        if (maxNx < Overlap.get(it, jt))
        {
            maxNx = it;
            maxNy = jt;
        }
    }

    if (maxNx == -1 or maxNy == -1)
    {
        maxOverlap.push_back(-1);
        maxOverlap.push_back(-1);
        maxOverlap.push_back(-1);
        return maxOverlap;
    }
    maxOverlap.push_back(maxNx);
    maxOverlap.push_back(maxNy);
    maxOverlap.push_back(Overlap.get(maxNx, maxNy));
    return maxOverlap;
}

Matrix<int> PhaseField::GetPhaseFieldOverlap(const int thPhase1,
        const int thPhase2)
{
    Matrix<int> Overlap;
    int numberOfGrains = FieldsStatistics.GrainStorage.size();
    Overlap.Allocate(numberOfGrains, numberOfGrains);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Fields,0,)
    {
        if (Interface(i,j,k))
        for (auto it = Fields(i,j,k).cbegin();
                  it < Fields(i,j,k).cend(); ++it)
        {
            int phaseIndex = it->index;
            int thPhaseIndex = FieldsStatistics[phaseIndex].Phase;
            if (thPhaseIndex == thPhase1)
            {
                for (auto jt = it+1;
                          jt < Fields(i,j,k).cend(); ++jt)
                {
                    int phaseIndex2 = jt->index;
                    int thPhaseIndex2 = FieldsStatistics[phaseIndex2].Phase;

                    #ifdef _OPENMP
                    #pragma omp critical
                    #endif
                    if (thPhaseIndex2 == thPhase2)
                    {
                        Overlap.add(phaseIndex, phaseIndex2, 1);
                    }
                }
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
    return Overlap;
}

void PhaseField::Advect(Velocities& Vel, BoundaryConditions& BC,
                                 double dt, int scheme)
{

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDot,0,)
    if (Fields(i,j,k).flag)
    {
        FieldsDot(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    switch(scheme)
    {
        case Upwind:
        {
            const double dx2 = 0.5/dx;
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDot,0,)
            if (Fields(i,j,k).flag)
            {
                int pInd = 0;
                int state = Liquid;
                for (auto alpha = Fields(i,j,k).begin();
                          alpha != Fields(i,j,k).end(); alpha++)
                if(FieldsStatistics[alpha->index].State > state)
                {
                    pInd = alpha->index;
                    state = FieldsStatistics[alpha->index].State;
                }

                pInd = FieldsStatistics[pInd].Phase;

                for (auto alpha = Fields(i,j,k).begin();
                          alpha != Fields(i,j,k).end(); alpha++)
                {
                    int ind = alpha->index;

                    FieldsDot(i,j,k).set(ind,
                     (dx2*((fabs(Vel.Average(i-1,j,k)[0]) + Vel.Average(i-1,j,k)[0])*Fields(i-1,j,k)[ind] +
                           (fabs(Vel.Average(i,j-1,k)[1]) + Vel.Average(i,j-1,k)[1])*Fields(i,j-1,k)[ind] +
                           (fabs(Vel.Average(i,j,k-1)[2]) + Vel.Average(i,j,k-1)[2])*Fields(i,j,k-1)[ind] +
                           (fabs(Vel.Average(i+1,j,k)[0]) - Vel.Average(i+1,j,k)[0])*Fields(i+1,j,k)[ind] +
                           (fabs(Vel.Average(i,j+1,k)[1]) - Vel.Average(i,j+1,k)[1])*Fields(i,j+1,k)[ind] +
                           (fabs(Vel.Average(i,j,k+1)[2]) - Vel.Average(i,j,k+1)[2])*Fields(i,j,k+1)[ind]) -
                           (fabs(Vel.Average(i,j,k)[0]) +
                            fabs(Vel.Average(i,j,k)[1]) +
                            fabs(Vel.Average(i,j,k)[2])) * Fields(i,j,k)[ind]/dx));
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            break;
        }
        case LaxWendroff:
        {
            const double dx4 = 0.25/dx;
            const double dtx = dt/dx;

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDot,0,)
            if (Fields(i,j,k).flag)
            {
                int ind = 0;
                int state = Liquid;

                for (auto alpha = Fields(i,j,k).begin();
                          alpha != Fields(i,j,k).end(); alpha++)
                if(FieldsStatistics[alpha->index].State > state)
                {
                    ind = alpha->index;
                    state = FieldsStatistics[alpha->index].State;
                }

                ind = FieldsStatistics[ind].Phase;

                for (auto alpha = Fields(i,j,k).begin();
                          alpha != Fields(i,j,k).end(); alpha++)
                {
                    double val   = Fields(i,j,k)[alpha->index];
                    double valXm = Fields(i-1,j,k)[alpha->index];
                    double valXp = Fields(i+1,j,k)[alpha->index];
                    double valYm = Fields(i,j-1,k)[alpha->index];
                    double valYp = Fields(i,j+1,k)[alpha->index];
                    double valZm = Fields(i,j,k-1)[alpha->index];
                    double valZp = Fields(i,j,k+1)[alpha->index];

                    double ux = Vel.Phase(i,j,k)({ind})[0];
                    double uy = Vel.Phase(i,j,k)({ind})[1];
                    double uz = Vel.Phase(i,j,k)({ind})[2];

                    double uxm = Vel.Phase(i-1,j,k)({ind})[0];
                    double uxp = Vel.Phase(i+1,j,k)({ind})[0];

                    double uym = Vel.Phase(i,j-1,k)({ind})[1];
                    double uyp = Vel.Phase(i,j+1,k)({ind})[1];

                    double uzm = Vel.Phase(i,j,k-1)({ind})[2];
                    double uzp = Vel.Phase(i,j,k+1)({ind})[2];

                    FieldsDot(i,j,k).set(alpha->index,
                    dx4*(
                    // Along x
                    (val + valXm - dtx*(val*ux - valXm * uxm))*
                    (ux + uxm) -
                    (val + valXp - dtx*(valXp * uxp - val*ux))*
                    (uxp + ux) +
                    // Along y
                    (val + valYm - dtx*(val*uy - valYm * uym))*
                    (uy + uym) -
                    (val + valYp - dtx*(valYp * uyp - val*uy))*
                    (uyp + uy) +
                    // Along z
                    (val + valZm - dtx*(val*uz - valZm * uzm))*
                    (uz + uzm) -
                    (val + valZp - dtx*(valZp * uzp - val*uz))*
                    (uzp + uz)));
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            break;
        }
        default:
        {
            Info::WriteExit("Wrong advection scheme selected",
                             thisclassname, "PhaseField::Advect(Vel, BC, dt, scheme)");
            exit(13);
        }
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,FieldsDot,0,)
    if (Fields(i,j,k).flag)
    {
        for (auto alpha = FieldsDot(i,j,k).begin();
                  alpha != FieldsDot(i,j,k).end(); alpha++)
        {
            Fields(i,j,k).add(alpha->index, alpha->value*dt);

        }
        FieldsDot(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Finalize(BC);
}


void PhaseField::Advect(Velocities& Vel, BoundaryConditions& BC, FlowSolverLBM& LBM,
                                 double dt, int scheme)
{
    LBM.DetectObstaclesSimple(*this);
    Advect(Vel, BC, dt, scheme);
    LBM.DetectObstacles(*this);
    LBM.SetFluidNodesNearObstacle();
}


void PhaseField::init_omp_locks()
{
#ifdef _OPENMP
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, Fields.Bcells(), )
    {
        omp_init_lock(&omp_locks(i, j, k));
    }
    OMP_PARALLEL_STORAGE_LOOP_END
#endif
}

void PhaseField::destroy_omp_locks()
{
#ifdef _OPENMP
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Fields, Fields.Bcells(), )
    {
        omp_destroy_lock(&omp_locks(i, j, k));
    }
    OMP_PARALLEL_STORAGE_LOOP_END
#endif
}
}// namespace openphase
