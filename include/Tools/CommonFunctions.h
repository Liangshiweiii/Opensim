#ifndef COMMONFUNCTIONS_H
#define COMMONFUNCTIONS_H

#include "Tools/ElasticityTensors.h"
namespace opensim
{
class CommonFunctions
{
 public:
     static void CalculateDistancePeriodic(dVector3 A, dVector3 B,
             dVector3 &dist, double Nx, double Ny, double Nz)
     {
     	dist[0] = A[0]-B[0];
     	if (std::fabs(dist[0]) > std::fabs(dist[0]+Nx)) dist[0] += Nx;
     	else if (std::fabs(dist[0]) > std::fabs(dist[0]-Nx)) dist[0] -= Nx;
     	dist[1] = A[1]-B[1];
     	if (std::fabs(dist[1]) > std::fabs(dist[1]+Ny)) dist[1] += Ny;
     	else if (std::fabs(dist[1]) > std::fabs(dist[1]-Ny)) dist[1] -= Ny;
     	dist[2] = A[2]-B[2];
     	if (std::fabs(dist[2]) > std::fabs(dist[2]+Nz)) dist[2] += Nz;
     	else if (std::fabs(dist[2]) > std::fabs(dist[2]-Nz)) dist[2] -= Nz;
     }
};

}// namespace opensim
#endif