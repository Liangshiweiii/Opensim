

#include "FluidDynamics/InteractionSolidFluid.h"
#include "PhaseField.h"
#include "GrainInfo.h"
#include "Velocities.h"
#include "Tools/EulerAngles.h"

using namespace std;
namespace opensim
{

void InteractionSolidFluid::CollectGrainsStatistics(PhaseField& Phase)
{
    CollectGrainsStatisticsStepOne(Phase);
    CollectGrainsStatisticsStepTwo(Phase);
}

void InteractionSolidFluid::CollectGrainsStatisticsStepOne(PhaseField& Phase)
{

    for(unsigned int idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        Phase.FieldsStatistics[idx].Rcm.set_to_zero();
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        if(Phase.Interface(i,j,k))
        {
            for(auto it = Phase.Fields(i,j,k).cbegin();
                     it < Phase.Fields(i,j,k).cend(); ++it)
            {
                #ifdef _OPENMP
                #pragma omp critical
                #endif
                {
                    Phase.FieldsStatistics[it->index].Rcm[0] += i*it->value;
                    Phase.FieldsStatistics[it->index].Rcm[1] += j*it->value;
                    Phase.FieldsStatistics[it->index].Rcm[2] += k*it->value;
                }
            }
        }
        else
        {
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                int locIndex = Phase.Fields(i,j,k).front().index;
                Phase.FieldsStatistics[locIndex].Rcm[0] += i;
                Phase.FieldsStatistics[locIndex].Rcm[1] += j;
                Phase.FieldsStatistics[locIndex].Rcm[2] += k;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    for(unsigned int idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        if(Phase.FieldsStatistics[idx].Volume > 0.0)
        {
            Phase.FieldsStatistics[idx].Rcm *= 1.0/double(Phase.FieldsStatistics[idx].Volume);
        }
    }
}

void InteractionSolidFluid::CalculateCenterOfMassWithPeriodicBoundaryConditions(PhaseField& Phase)
{
    int Nx = Phase.Nx;
    int Ny = Phase.Ny;
    int Nz = Phase.Nz;
    GrainInfo gloFieldsStatistics1;
    GrainInfo gloFieldsStatistics2;
    gloFieldsStatistics1.Allocate(Phase.FieldsStatistics.size());
    gloFieldsStatistics2.Allocate(Phase.FieldsStatistics.size());
    for(unsigned int idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        gloFieldsStatistics1[idx].Rcm.set_to_zero();
        gloFieldsStatistics2[idx].Rcm.set_to_zero();
    }
    //#pragma omp parallel //OMP BEGIN
    {
        GrainInfo locFieldsStatistics1;
        GrainInfo locFieldsStatistics2;
        locFieldsStatistics1.Allocate(Phase.FieldsStatistics.size());
        locFieldsStatistics2.Allocate(Phase.FieldsStatistics.size());

        //#pragma omp for collapse(OMP_COLLAPSE_LOOPS) schedule(dynamic,OMP_DYNAMIC_CHUNKSIZE)
        STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
        {
            if(Phase.Interface(i,j,k))
            {
                for(auto it = Phase.Fields(i,j,k).cbegin();
                         it < Phase.Fields(i,j,k).cend(); ++it)
                {
                    double theta = (double)i/(double)Nx*2.0*Pi;
                    locFieldsStatistics1[it->index].Rcm[0] += Nx/(2.0*Pi)*cos(theta)*it->value;
                    locFieldsStatistics2[it->index].Rcm[0] += Nx/(2.0*Pi)*sin(theta)*it->value;
                    theta = (double)j/(double)Ny*2.0*Pi;
                    locFieldsStatistics1[it->index].Rcm[1] += Ny/(2.0*Pi)*cos(theta)*it->value;
                    locFieldsStatistics2[it->index].Rcm[1] += Ny/(2.0*Pi)*sin(theta)*it->value;
                    theta = (double)k/(double)Nz*2.0*Pi;
                    locFieldsStatistics1[it->index].Rcm[2] += Nz/(2.0*Pi)*cos(theta)*it->value;
                    locFieldsStatistics2[it->index].Rcm[2] += Nz/(2.0*Pi)*sin(theta)*it->value;
                }
            }
            else
            {
                int locIndex = Phase.Fields(i,j,k).front().index;
                double theta = (double)i/(double)Nx*2.0*Pi;
                locFieldsStatistics1[locIndex].Rcm[0] += Nx/(2.0*Pi)*cos(theta);
                locFieldsStatistics2[locIndex].Rcm[0] += Nx/(2.0*Pi)*sin(theta);
                theta = (double)j/(double)Ny*2.0*Pi;
                locFieldsStatistics1[locIndex].Rcm[1] += Ny/(2.0*Pi)*cos(theta);
                locFieldsStatistics2[locIndex].Rcm[1] += Ny/(2.*Pi)*sin(theta);
                theta = (double)k/(double)Nz*2.0*Pi;
                locFieldsStatistics1[locIndex].Rcm[2] += Nz/(2.0*Pi)*cos(theta);
                locFieldsStatistics2[locIndex].Rcm[2] += Nz/(2.0*Pi)*sin(theta);
            }
        }
        STORAGE_LOOP_END
        //#pragma omp critical
        {
            for(unsigned int idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
            {
                gloFieldsStatistics1[idx].Rcm += locFieldsStatistics1[idx].Rcm;
                gloFieldsStatistics2[idx].Rcm += locFieldsStatistics2[idx].Rcm;
            }
            for(unsigned int idx = 0; idx < Phase.FieldsStatistics.size(); idx++) {
                if(Phase.FieldsStatistics[idx].Volume > 0.0)
                {
                    gloFieldsStatistics1[idx].Rcm *= 1.0/double(Phase.FieldsStatistics[idx].Volume);
                    gloFieldsStatistics2[idx].Rcm *= 1.0/double(Phase.FieldsStatistics[idx].Volume);
                }
                Phase.FieldsStatistics[idx].Rcm[0] = atan2(-gloFieldsStatistics2[idx].Rcm[0],-gloFieldsStatistics1[idx].Rcm[0]) + Pi;
                Phase.FieldsStatistics[idx].Rcm[0] = Nx*Phase.FieldsStatistics[idx].Rcm[0]/(2.0*Pi);
                Phase.FieldsStatistics[idx].Rcm[1] = atan2(-gloFieldsStatistics2[idx].Rcm[1],-gloFieldsStatistics1[idx].Rcm[1]) + Pi;
                Phase.FieldsStatistics[idx].Rcm[1] = Nx*Phase.FieldsStatistics[idx].Rcm[1]/(2.0*Pi);
                Phase.FieldsStatistics[idx].Rcm[2] = atan2(-gloFieldsStatistics2[idx].Rcm[2],-gloFieldsStatistics1[idx].Rcm[2]) + Pi;
                Phase.FieldsStatistics[idx].Rcm[2] = Nx*Phase.FieldsStatistics[idx].Rcm[2]/(2.0*Pi);
            }
        }
    }// OMP END

}

void InteractionSolidFluid::CollectGrainsStatisticsStepTwo(PhaseField& Phase)
{
    for(unsigned int idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        Phase.FieldsStatistics[idx].InertiaM.set_to_zero();
    }
    //#pragma omp parallel
    {
        GrainInfo locFieldsStatistics;
        locFieldsStatistics.Allocate(Phase.FieldsStatistics.size());

        //#pragma omp for collapse(OMP_COLLAPSE_LOOPS) schedule(dynamic,OMP_DYNAMIC_CHUNKSIZE)
        STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
        {
            if(Phase.Interface(i,j,k))
            {
                for(auto it = Phase.Fields(i,j,k).cbegin();
                         it != Phase.Fields(i,j,k).cend(); ++it)
                {
                    double locRadX = (i - Phase.FieldsStatistics[it->index].Rcm[0])*Phase.dx;
                    double locRadY = (j - Phase.FieldsStatistics[it->index].Rcm[1])*Phase.dx;
                    double locRadZ = (k - Phase.FieldsStatistics[it->index].Rcm[2])*Phase.dx;

                    locFieldsStatistics[it->index].InertiaM(0,0) += it->value*(locRadY*locRadY + locRadZ*locRadZ);
                    locFieldsStatistics[it->index].InertiaM(1,1) += it->value*(locRadX*locRadX + locRadZ*locRadZ);
                    locFieldsStatistics[it->index].InertiaM(2,2) += it->value*(locRadX*locRadX + locRadY*locRadY);

                    locFieldsStatistics[it->index].InertiaM(0,1) -= it->value*(locRadX*locRadY);
                    locFieldsStatistics[it->index].InertiaM(1,0) -= it->value*(locRadY*locRadX);
                    locFieldsStatistics[it->index].InertiaM(0,2) -= it->value*(locRadX*locRadZ);
                    locFieldsStatistics[it->index].InertiaM(2,0) -= it->value*(locRadZ*locRadX);
                    locFieldsStatistics[it->index].InertiaM(1,2) -= it->value*(locRadY*locRadZ);
                    locFieldsStatistics[it->index].InertiaM(2,1) -= it->value*(locRadZ*locRadY);
                }
            }
            else
            {
                int locIndex = Phase.Fields(i, j, k).front().index;
                double locRadX = (i - Phase.FieldsStatistics[locIndex].Rcm[0])*Phase.dx;
                double locRadY = (j - Phase.FieldsStatistics[locIndex].Rcm[1])*Phase.dx;
                double locRadZ = (k - Phase.FieldsStatistics[locIndex].Rcm[2])*Phase.dx;

                locFieldsStatistics[locIndex].InertiaM(0,0) += (locRadY*locRadY + locRadZ*locRadZ);
                locFieldsStatistics[locIndex].InertiaM(1,1) += (locRadX*locRadX + locRadZ*locRadZ);
                locFieldsStatistics[locIndex].InertiaM(2,2) += (locRadX*locRadX + locRadY*locRadY);

                locFieldsStatistics[locIndex].InertiaM(0,1) -= (locRadX*locRadY);
                locFieldsStatistics[locIndex].InertiaM(1,0) -= (locRadY*locRadX);
                locFieldsStatistics[locIndex].InertiaM(0,2) -= (locRadX*locRadZ);
                locFieldsStatistics[locIndex].InertiaM(2,0) -= (locRadZ*locRadX);
                locFieldsStatistics[locIndex].InertiaM(1,2) -= (locRadY*locRadZ);
                locFieldsStatistics[locIndex].InertiaM(2,1) -= (locRadZ*locRadY);
            }
        }
        STORAGE_LOOP_END
        //#pragma omp critical
        {
            for(unsigned int idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
            {
                Phase.FieldsStatistics[idx].InertiaM += locFieldsStatistics[idx].InertiaM;
            }
        }
    }// OMP END

    double dx = Phase.dx;
    double dV = dx*dx*dx;
    for(unsigned int idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    {
        Phase.FieldsStatistics[idx].InertiaM *= dV*Phase.FieldsStatistics[idx].Density;
    }
}

void InteractionSolidFluid::CalculateSolidVelocities(PhaseField& Phase, Velocities& Vel, double dt)
{
    double dx = Phase.dx;
    double dV = dx*dx*dx;
    CollectGrainsStatistics(Phase);
    for(size_t idx = 0; idx < Phase.FieldsStatistics.size(); idx++)
    if(Phase.FieldsStatistics[idx].Exist &&
       Phase.FieldsStatistics[idx].State == Solid &&
       Phase.FieldsStatistics[idx].IsMobile)
    {
        double SolidMass_1 = 1.0/(Phase.FieldsStatistics[idx].Volume*dV*Phase.FieldsStatistics[idx].Density);
        Phase.FieldsStatistics[idx].Acm  = Phase.FieldsStatistics[idx].Force * SolidMass_1;
        Phase.FieldsStatistics[idx].Vcm += Phase.FieldsStatistics[idx].Acm * dt;
        Phase.FieldsStatistics[idx].aAcc = Phase.FieldsStatistics[idx].InertiaM.inverted() * Phase.FieldsStatistics[idx].Torque;
        Phase.FieldsStatistics[idx].aVel += Phase.FieldsStatistics[idx].aAcc * dt;
        EulerAngles locAngles({Phase.FieldsStatistics[idx].aVel[0] * dt,
                               Phase.FieldsStatistics[idx].aVel[1] * dt,
                               Phase.FieldsStatistics[idx].aVel[2] * dt}, XYZ);
        Phase.FieldsStatistics[idx].Orientation += locAngles.getQuaternion();
    }

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0,)
    {
        for(auto alpha = Phase.Fields(i,j,k).cbegin();
                 alpha != Phase.Fields(i,j,k).cend(); ++alpha)
        if (Phase.FieldsStatistics[alpha->index].State == Solid &&
            Phase.FieldsStatistics[alpha->index].IsMobile)
        {
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

            dVector3 R;
            R[0] = (i - Phase.FieldsStatistics[alpha->index].Rcm[0])*dx;
            R[1] = (j - Phase.FieldsStatistics[alpha->index].Rcm[1])*dx;
            R[2] = (k - Phase.FieldsStatistics[alpha->index].Rcm[2])*dx;

            Vel.Phase(i,j,k)({Phase.FieldsStatistics[alpha->index].Phase}) = Phase.FieldsStatistics[alpha->index].Vcm + W*R;
        }
        else if (Phase.FieldsStatistics[alpha->index].State == Solid)
        {
            Vel.Phase(i,j,k)({Phase.FieldsStatistics[alpha->index].Phase}) = Phase.FieldsStatistics[alpha->index].Vcm;
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END

    Vel.CalculateAverage(Phase);
}

} //namespace opensim
