#ifndef INTERFACEENERGY_H
#define INTERFACEENERGY_H

#include "Tools/Includes.h"
#include "Tools/Node.h"
#include "Info.h"
#include "Settings.h"

namespace opensim{
    class ElasticProperties;
    class PhaseField;

    class InterfaceEnergy : public OPobject{
        public:
        int Nphases;
        int Nx;
        int Ny;
        int Nz;
        double iWidth;
        bool Averaging;

        double maxSigma;

        Matrix< double > IntEnergy;
        Matrix< double > MinIntEnergy;
        Matrix< double > Anisotropy;
        Tensor< double, 3 > Tau;
        Storage3D< Node, 0 > IntEnergy3D;
        Storage3D< Node, 0 > RawIntEnergy3D;
        Storage3D< Node, 0 > Counter;

        InterfaceEnergy(const Settings& locSettings, unsigned int boundary = 1);
        InterfaceEnergy(void)
        {

        }
        using OPObject::ReadInput;
        void ReadInput(std::string InputFileName);

        int Bcells(void) const
        {
            return IntEnergy3D.Bcells();
        }

        double operator()(const int i, const int j, const int k, const int alpha, const int beta) const
        {
            return IntEnergy3D(i,j,k).get_sym(alpha, beta);
        };
        void set(const int i, const int j, const int k, const int alpha, const int beta, const double value)
        {
            IntEnergy3D(i,j,k).set_sym(alpha, beta, value);
        };
        void add(const int i, const int j, const int k, const int alpha, const int beta, const double value)
        {
            IntEnergy3D(i,j,k).add_sym(alpha, beta, value);
        };
        void clear(const int i, const int j, const int k)
        {
            IntEnergy3D(i,j,k).clear();
            RawIntEnergy3D(i,j,k).clear();
            Counter(i,j,k).clear();
        };
        void set_raw(const int i, const int j, const int k, const int alpha, const int beta, const double value)
        {
            RawIntEnergy3D(i,j,k).set_sym(alpha, beta, value);
        };
        void add_raw(const int i, const int j, const int k, const int alpha, const int beta, const double value)
        {
            RawIntEnergy3D(i,j,k).add_sym(alpha, beta, value);
        };
        void Initialize(const Settings& locSettings, unsigned int boundary = 1);
        void ReInitialize(const PhaseField& Phase);
        void CalculateCubic(const PhaseField& Phase);
        void CalculateHex(const PhaseField& Phase);
        void ConsiderStrain(const PhaseField& Phase, const ElasticProperties& EP);
        void Set(const PhaseField& Phase);
        void SetRaw(const PhaseField& Phase);
        void Average(const PhaseField& Phase, const BoundaryConditions& BC);

        /**************************************************************************/
        void WriteVTK(const int tStep, const int indexA, const int indexB);
        void WriteVTK(const int tStep, const PhaseField& Phase, const int indexA);
        /**************************************************************************/
    };
}
#endif