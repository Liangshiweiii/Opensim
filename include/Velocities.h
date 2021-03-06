#ifndef VELOCITIES_H
#define VELOCITIES_H

#include "Tools/Includes.h"
namespace opensim
{

class Settings;
class PhaseField;
class BoundaryConditions;

class Velocities : public OPObject                                              ///< Storage for the advection velocities
{
 public:
    Velocities(){};
    Velocities(const Settings& locSettings, const unsigned int boundary = 2);

    void Initialize(const Settings& locSettings,
            const unsigned int boundary = 2);                                   ///<  Allocates memory, initializes the settings
    void SetBoundaryConditions(const BoundaryConditions& BC);                   ///<  Sets boundary conditions for the velocities
    void SetAverage(dVector3& value);                                           ///<  Set average velocity to a fixed value
    void SetAllPhases(dVector3& value);                                         ///<  Set phase velocities to a fixed value
    void CalculateAverage(const PhaseField& Phase);                             ///<  Calculates the average velocity in each grid point using phase velocities
    void PrescribePhaseVelocities(PhaseField& Phase);                           ///<  Sets phase velocities equal to average velocities
    void Clear(void);                                                           ///<  Clears the velocity storage
    void WriteAverageVTK(int tStep) const;                                      ///<  Writes average velocities in the VTK format into a file
    void WritePhaseVTK(int tStep, int idx) const;                               ///<  Writes phase velocities in the VTK format into a file
    void WriteStatistics(int tStep, double dt);                                 ///<  Writes velocity statistics into a file, input: time step
    double getMaxVelocity();                                                    ///<  returns maximal velocity component
    double getMaxVelocityOMP();                                                 ///<  returns maximal velocity component
    dVector3 getCriticalPoint();                                                ///<  returns the coordinates of the point with the highest velocity
    void Remesh(int newNx, int newNy, int newNz, BoundaryConditions& BC);       ///<  Remesh the system while keeping the data
    void MoveFrame(int dx, int dy, int dz, BoundaryConditions& BC);             ///<  Moves frame by (dx,dy,dz).
    void Advect(BoundaryConditions& BC, double dt, int scheme = Upwind);        ///<  Advects average velocities
    
    void PrintPointStatistics(int x, int y, int z);                             ///<  Prints local velocity values to screen.
    Storage3D < dVector3, 1 > Phase;                                            ///<  Phase velocities storage
    Storage3D < dVector3, 0 > Average;                                          ///<  Average velocities storage
    Storage3D < dVector3, 0 > AverageDot;                                       ///<  Average velocities increments storage
    int Nx;                                                                     ///<  Size of the inner calculation domain along X
    int Ny;                                                                     ///<  Size of the inner calculation domain along Y
    int Nz;                                                                     ///<  Size of the inner calculation domain along Z
    double dx;                                                                  ///<  Grid spacing
    int Nphases;                                                                ///<  Number of thermodynamic phases

    Velocities& operator= (const Velocities& rhs);                              ///< Copy operator for Velocities class

 protected:
 private:
    dVector3 criticalPoint;                                                     ///<  Stores the coordinates of the point with the highest velocity
    double   maxVelocity;
};

}
#endif