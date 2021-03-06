#ifndef INTERFACEMOBILITY_H
#define INTERFACEMOBILITY_H

#include "Tools/Node.h"
#include "Tools/OPObjects.h"
#include "Info.h"
#include "Settings.h"

namespace opensim
{
class BoundaryConditions;
class PhaseField;
class Temperature;

class InterfaceMobility : public OPObject                                       /// 界面移动性模块的纯虚拟“接口”。Pure virtual "interface" for the Interfacial mobility module.
{
 public:
    InterfaceMobility(){};
    InterfaceMobility(const Settings& locSettings, unsigned int boundary = 1);
    void Initialize(const Settings& locSettings, unsigned int boundary = 1);
    using OPObject::ReadInput;
    void ReadInput(const std::string FileName);                                 ///< 读取内部参数的初始值。Reads initial values for internal parameters.
    
    double operator()(const int i, const int j, const int k, const int alpha, const int beta)
    {
        return IntMobility3D(i,j,k).get_sym(alpha, beta);
    };
    void set(const int i, const int j, const int k, const int alpha, const int beta, const double value)
    {
        IntMobility3D(i,j,k).set_sym(alpha, beta, value);
    };
    void set_raw(const int i, const int j, const int k, const int alpha, const int beta, const double value)
    {
        RawIntMobility3D(i,j,k).set_sym(alpha, beta, value);
    };
    void add(const int i, const int j, const int k, const int alpha, const int beta, const double value)
    {
        IntMobility3D(i,j,k).add_sym(alpha, beta, value);
    };
    void add_raw(const int i, const int j, const int k, const int alpha, const int beta, const double value)
    {
        RawIntMobility3D(i,j,k).add_sym(alpha, beta, value);
    };
    void clear(const int i, const int j, const int k)
    {
        IntMobility3D(i,j,k).clear();
        RawIntMobility3D(i,j,k).clear();
        Counter(i,j,k).clear();
    };

    void ReInitialize(const PhaseField& Phase);
    void Set(const PhaseField& Phase);
    void Set(const PhaseField& Phase, const Temperature& Tx);
    void CalculateCubic(const PhaseField& Phase);
    void CalculateHex(const PhaseField& Phase);
    void CalculateCubic(const PhaseField& Phase, const Temperature& Tx);
    void CalculateHex(const PhaseField& Phase, const Temperature& Tx);
    void Average(const PhaseField& Phase, const BoundaryConditions& BC);
    void WriteVTK(const int tStep,const int indexA, const int indexB);    
    
    int Nx;
    int Ny;
    int Nz;
    double iWidth;
    double R;

    double iMu;
    int Nphases;
    bool Averaging;
    bool GrainBoundariesAreMobile;                                              /// 用于打开/关闭给定相的晶界迁移Used to turn on/off grain boundary migration for a given phase
    double maxMu;

    Matrix< double > IntMobility;
    Matrix< double > IntAnisotropy;
    Matrix< double > MaxIntMobility;
    Matrix< double > ActivationEnergy;
    Storage3D< Node, 0 > IntMobility3D;
    Storage3D< Node, 0 > RawIntMobility3D;
    Storage3D< Node, 0 > Counter;
};

}//namespace opensim
#endif  // header guard
