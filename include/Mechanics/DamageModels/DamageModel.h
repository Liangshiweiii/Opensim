

#ifndef DAMAGEMODEL_H
#define DAMAGEMODEL_H

#include "fftw3.h"
#include "Tools/Includes.h"

namespace opensim
{
class Settings;
class PhaseField;
class ElasticProperties;
class DamageProperties;

class DamageModel : public OPObject
{
 public:
    using OPObject::ReadInput;
    void ReadInput(const std::string IFileName);

    DamageModel(){};
    DamageModel(const Settings& locSettings, const std::string InputFileName = "NONE");
    ~DamageModel(void);

    void Initialize(const Settings& locSettings);
    void Reinitialize(const long int newNx, const long int newNy,
                                                          const long int newNz);

    void Solve(ElasticProperties& EP, DamageProperties& DP, PhaseField& Phi, BoundaryConditions& BC);

    void ExecuteForwardTRAFO_KAPPA(void);
    void ExecuteBackwardTRAFO_KAPPANL(void);
    void F_CalculateNLEquationPerlings(ElasticProperties& EP);
    void F_EValuateDamageSurfacePerlingsDuctile2Intf(PhaseField& Phi, ElasticProperties& EP,DamageProperties& DP);    ///correctly accounts for bulk and interface regions
    void F_NormalizeKAPPANL();                                                  /// Normalizes EEQNL after backward transformation
    void F_CombineBrittleDuctileDamage(DamageProperties& DP);								/// Adds up brittle and ductile damage to effectife damage
    void F_SetEquivalentPlasticStrain(DamageProperties& DP);
    void WriteVTK(int tStep);

    fftw_plan  ForwardPlanKAPPA;                                                ///
    fftw_plan  BackwardPlanKAPPANL;                                             ///

    double                 maxDmg;
    double               * damage;                                              /// Damage parameter Peerlings
    double               * kappa;
    double               * kappaNL;
    std::complex<double> * kappaF;
    std::complex<double> * kappaNLF;

 protected:

 private:
    int Nx;                                                                     /// System size along X direction
    int Ny;                                                                     /// System size along X direction
    int Nz;                                                                     /// System size along Z direction
    int Nzc;                                                                    /// System size along Fourier Z direction
    int rlSIZE;                                                                 /// System size (real space): nX*nY*nZ
    int rcSIZE;                                                                 /// System size (Fourier space): nX*nY*(nZ/2+1)
    int Nphases;
    double Norm;                                                                /// Normalization Coefficient: 1.0/rlSIZE
    double iLx;                                                                 /// 1.0/nX*dx constant
    double iLy;                                                                 /// 1.0/nY*dx constant
    double iLz;                                                                 /// 1.0/nZ*dx constant
    double dx;

    double * Q[3];                                                              /// Wave vectors

    Storage<bool> damageflag;
    double  alphaNL;
    Storage<double> k0_linear;
    Storage<double> kc_linear;

    void QXYZ(void);                                                            /// Sets the wave vectors

};
}// namespace opensim

#endif
