#ifndef DOUBLEOBSTACLE_H
#define DOUBLEOBSTACLE_H

#include "Tools/Includes.h"

namespace opensim{
    class BoundaryConditions;
    class InterfaceEnergy;
    class InterfaceMobility;
    class PhaseField;
    class Settings;

    class DoubleObstacle : public OPObject{
        public:
        DoubleObstacle(){};
        DoubleObstacle(Settings& locSettings);                                      ///< 用构造函数初始化模块Initializes module with constructor

        void Initialize(Settings& locSettings);                                     ///< 初始化，仅指示已创建模块。Initializes, just to indicate that module has been created.
        void CalculatePhaseFieldIncrements(PhaseField& Phase, InterfaceEnergy& Sigma,
                                                            InterfaceMobility& Mu);///< 计算与界面曲率有关的驱动力。Calculates interface curvature related driving force.
        void FixSpreading(PhaseField& Phase, BoundaryConditions& BC, double cutoff);///< 移除低于临界值的相场 (用于消除平流中的寄生扩散)
        double Energy(PhaseField& Phase, InterfaceEnergy& Sigma);                   ///< 提供模拟域中的总界面能。Provides total interfacial energy in the simulation domain.
        double AverageEnergyDensity(PhaseField& Phase, InterfaceEnergy& Sigma);     ///< 提供模拟域中的平均界面能量密度。Provides the average interfacial energy density in the simulation domain.
        double PointEnergy (PhaseField& Phase, InterfaceEnergy& Sigma,
                int i, int j, int k) ;                                              ///< 提供给定点（i，j，k）的界面能。Provides the interfacial energy in a given point (i,j,k).
        void WriteEnergyVTK(PhaseField& Phase, InterfaceEnergy& Sigma, int tStep);  ///< 在给定的时间步长tStep中以VTK格式写入界面能量Writes interfacial energy in VTK format for a given timestep tStep

        
    };
}

#endif