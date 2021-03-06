#ifndef TOOLS_H
#define TOOLS_H

#include "Tools/Includes.h"
#include "Tools/ElasticityTensors.h"
namespace opensim
{

class PhaseField;
class Orientations;

class Tools
{
 public:
    /*static void WriteAbaqusInput(const PhaseField& Phi,                    /// Writes Simulia Abaqus input file (Maintainer: PE)
                                    const Orientations& OR);
    */
    static void Decompose(dMatrix3x3& Deformation, dMatrix3x3& Rot,
                                                            dMatrix3x3& Strain);
    static void Eigensystem(dMatrix3x3& M, dMatrix3x3& Eigenvectors,
                                                       dMatrix3x3& Eigenvalues);// M MUST BE _SYMMETRIC_, calculates eigenvalues and (column) eigenvectors of M
    static void Eigensystem2(dMatrix3x3& M, dMatrix3x3& Eigenvectors,
                                                       dMatrix3x3& Eigenvalues, bool calculateEigenVectors = true); // Non-iterative! M MUST BE _SYMMETRIC_, calculates eigenvalues and (column) eigenvectors of M
    /**************************************************************************/

 protected:
 private:
};
}
#endif