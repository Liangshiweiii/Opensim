#ifndef ORIENTATIONS_H
#define ORIENTATIONS_H

#include "Tools/Includes.h"
#include "PhaseField.h"

namespace opensim
{

class Settings;
class BoundaryConditions;
class PhaseField;
class Quaternion;
class Crystallography;
class Velocities;

class Orientations : public OPObject                                            ///< Module which stores and handles orientation data like rotation matrices, Euler angles etc.
{

 public:

    Orientations(){};
    Orientations(Settings& locSettings);

    void Initialize(const Settings& locSettings);                               ///< Initializes the module, allocate the storage, assign internal variables

    void SetBoundaryConditions(const BoundaryConditions& BC);                   ///< Sets boundary conditions
    void SetRandomGrainOrientations(PhaseField& Phase, const int seed = 0);     ///< Sets random grains orinetations
    void Remesh(const int newNx, const int newNy, const int newNz,
                                                        BoundaryConditions& BC);///< Remesh and reallocate orientations
    void WriteVTK(const int tStep);                                             ///< Writes rotations to the VTK file
    void WriteTotalVTK(PhaseField& Phase, const int tStep) const;               ///< Writes total rotations (including grains orientations) to the VTK file
    void Write(const int tStep) const;                                          ///< Writes orientations storage to binary file
    void Read(const int tStep);                                                 ///< Reads orientations from binary file
    void WriteGrainEBSDDataQuaternions(const PhaseField& Phase, const int tStep);///< Writes EBSD file (for MTex). In interface, majority phase is used.
    void WriteMisorientationsVTK(const int tStep, const std::string measure) const;
    void WriteRotated100Vector(const int tStep);
    static double getMisorientation(const dMatrix3x3 RotMatA, const dMatrix3x3 RotMatB); ///< Calculates missorientation between two matrices without consideration of symmetries
    static double getMisorientationCubic(const dMatrix3x3 RotMatA, const dMatrix3x3 RotMatB, const Crystallography& CR); /// Calculates missorientation between two matrices without consideration of symmetries

    void PrintPointStatistics(int x, int y, int z);
    void Advect(Velocities& Vel, BoundaryConditions& BC, double dt, int scheme = Upwind); ///< Advects orientations
    int Nx;
    int Ny;
    int Nz;

    double dx;

    int Nphases;

    Storage3D<Quaternion, 0>    Quaternions;
    Storage3D<Quaternion, 0>    QuaternionsDot;

    Orientations& operator= (const Orientations& rhs);

    dMatrix3x3 getEffectiveGrainRotation(const PhaseField& Phase,
                                    const int i, const int j, const int k) const///< Returns averaged grain rotations, averaging is done via quaternions
    {
        Quaternion locQuaternion;

        if(Phase.Fields(i,j,k).flag < 2)
        {
            int index = Phase.Fields(i,j,k).cbegin()->index;
            locQuaternion = Phase.FieldsStatistics[index].Orientation;
        }
        else
        {
            locQuaternion.set(0,0,0,0);

            for(auto alpha = Phase.Fields(i,j,k).cbegin();
                     alpha < Phase.Fields(i,j,k).cend(); alpha++)
            {
                locQuaternion += Phase.FieldsStatistics[alpha->index].Orientation*alpha->value;
            }
        }

        return locQuaternion.RotationMatrix;
    }

    dMatrix3x3 getTotalRotation(const PhaseField& Phase,
                                    const int i, const int j, const int k) const///< Returns the local rotation including grain orientation and deformation induced rotation
    {
        return Quaternions(i,j,k).RotationMatrix*getEffectiveGrainRotation(Phase, i,j,k);
    }

    dMatrix3x3 getTotalGrainRotation(PhaseField& Phase, int i, int j, int k,
                                                                      int alpha)///< Returns the total rotation seen by the specified grain including deformation induced rotation
    {
        return Quaternions(i,j,k).RotationMatrix*Phase.FieldsStatistics[alpha].Orientation.RotationMatrix;
    }

 protected:
 private:
};
}// namespace opensim
#endif
