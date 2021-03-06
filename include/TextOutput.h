#ifndef TEXTOUTPUT_H
#define TEXTOUTPUT_H

#include "Tools/Includes.h"

namespace opensim
{
class BoundaryConditions;
class ChemicalProperties;
class ElasticProperties;
class GrainInfo;
class InterfaceField;
class PhaseField;
class PhaseField;
class PlasticFlowCP;
class Settings;
class Temperature;
class Composition;
class ThermodynamicProperties;

class TextOutput : public OPObject                                              ///< 列表模拟数据的输出
{
 public:
    static void PhasePercent(ChemicalProperties& CP, PhaseField& Phi,
                               Settings& OPSettings, std::string filename,
                               double time);                                    ///< 写入热力学相的体积百分比
    static void GrainVolumes(PhaseField& Phi, std::string filename, double time);    ///< 写入每个晶粒的体积
    static void LineConcentration(Composition& Cx, ChemicalProperties& CP, PhaseField& Phi,
                                  std::string filename, double timestep,
                                  std::string type, std::string axis, int x, int y,int z);///< 在分离的文件中将成分写成一条直线
    static void AverageStress(ElasticProperties& EP, std::string filename,
                              double timeOrStrain);                             ///< 记录随时间/应变的平均应力
    static void AverageStrain(Storage3D<vStrain,0>& Strains, std::string filename,
                              double timeOrStrain);                             ///< 记录一段时间内的平均应变
    static void AverageDouble(Storage3D<double,0>& value, std::string filename,
                              double timeOrStrain);
    static void AverageCRSS(PlasticFlowCP& PFCP, PhaseField& Phase,
                            std::string filename, double timeOrStrain);         ///< 编写随时间/应变的平均CRSS
    static void AverageTemp(Temperature& Tx, std::string filename, double time);///< 写入平均系统温度
 private:
    static bool FileExists(std::string filename);                               ///< 检查给定文件名的可用性
    static void LocalPhaseComposition(Composition& Cx, ChemicalProperties& CP,
                                      PhaseField& Phi, std::string filename,
                                      double time, int x, int y, int z);         ///<
};

} // namespace opensim
#endif
