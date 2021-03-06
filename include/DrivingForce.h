#ifndef DRIVINGFORCE_H
#define DRIVINGFORCE_H

#include "Tools/Includes.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"

namespace opensim
{
    class Settings;
    class Node;
    class PhaseField;
    class InterfaceMobility;
    class InterfaceEnergy;
    class BoundaryConditions;

    class DrivingForce : public OPObject
    {
        public:
        DrivingForce() : Averaging(0) {};
        DrivingForce(const Settings& locSettings, const int boundary = 1);          ///< 初始化驱动力类别的存储和内部变量。Initializes the storage and internal variables of the driving force class.
        void Initialize(const Settings& locSettings, const int boundary = 1);       ///< 初始化驱动力类别的存储和内部变量。Initializes the storage and internal variables of the driving force class.
        using OPObject::ReadInput;
        void ReadInput(const std::string InputFileName);                            ///< 读取驱动力设置Reads driving force settings

        void Clear();                                                               ///< 清空存储空间Empties the storage
        void Average(PhaseField& Phase, BoundaryConditions& BC);                    ///< 界面上的平均驱动力
        void CollectAverage(PhaseField& Phase);                                     ///< 平均值的第一部分First part of Average
        void DistributeAverage(PhaseField& Phase);                                  ///< 平均值的第二部分Second part of Average
        void SetBoundaryConditions(BoundaryConditions& BC);                         ///< 设定边界条件Sets the boundary conditions
        void MergePhaseFieldIncrements(PhaseField& Phase, InterfaceEnergy& Sigma,
                InterfaceMobility& Mu);                                             ///< 将驱动力合并到相场增量中。Merges the driving force into the phase field increments.
        double MaxTimeStep(PhaseField& Phase, InterfaceEnergy& Sigma,
                InterfaceMobility& Mu, Settings& OPSettings, double TheorLimit,
                double NumLimit);                                                   ///< 计算给定时间步的最小相场迭代次数Calculates minimum number of phase field iterations for a given time step
        void AddNoiseRelative(double magnitude, int alpha, int beta);               ///< 为alpha和beta对添加噪声与实际驱动力的相对大小adding noise with a relative magnitude of the actual driving force for the pair of alpha and beta
        void AddNoiseAbsolute(double amplidude, int alpha, int beta);               ///< 添加与本地驱动力无关的噪声，由绝对振幅选择adding noise independent of local driving force, chosen by absolute amplitude
        void PrintDiagnostics(void);                                                ///< 将驱动力统计信息打印到屏幕上。指示推力是否超调。Prints the driving force statistics to the screen. Indicate if the riving force is overshooting.
        void PrintPointStatistics(const int i, const int j, const int k) const;
        int Nx;                                                                     ///< 系统的X维度。X dimension of the system.
        int Ny;                                                                     ///< 系统的Y维度。Y dimension of the system.
        int Nz;                                                                     ///< 系统的Z维度。Z dimension of the system.
        int Range;                                                                  ///< 给定网格点周围球体的半径，在该球体上平均驱动力Radius of a sphere around the given grid point over which the driving force is averaged
        double PhiThreshold;                                                        ///< 概述用于驱动力平均的界面的内部。Outlines the inner part of the interface for driving force averaging.
        void WriteVTK(const int tStep, const int indexA, const int indexB) const;   ///< 将从具有索引B的相场作用到具有索引A的相场的驱动力Writes the driving force acting from phase field with indexB onto phase field with indexA
        void WriteVTKforPhases(PhaseField& Phi, const int tStep) const;             ///< 写出作用在所有热力学相之间而不是单个晶粒之间的平均驱动力Writes the average driving force acting between all thermodynamic phases, not between individual grains
        void Remesh(int newNx, int newNy, int newNz);                               ///< 在保留数据的同时重新划分存储空间Remesh the storage while keeping the data
        Storage3D<Node, 0> Raw;                                                     ///< 原始驱动力存储Raw driving force storage
        //experimental
        std::mt19937_64 NoiseGenerator;
        std::normal_distribution <double> NoiseDistribution;
        //end experimental
        DrivingForce& operator= (const DrivingForce& rhs);                          ///< DrivingForce类的复制运算符Copy operator for DrivingForce class
        
        protected:
        int dGabLimitReachedCounter;                                                ///< 累积两次调用PrintDiagnostics（）之间的驱动力超调事件Accumulates the driving force overshooting events between calls to PrintDiagnostics()
        double MAXdGabOvershoot;                                                    ///< 保持在两次调用PrintDiagnostics（）之间观察到的最大驱动力超调Holds the maximum driving force overshoot observed between calls to PrintDiagnostics()
        double MAXDrivingForce;                                                     ///< 保持在两次调用PrintDiagnostics（）之间观察到的最大驱动力Holds the maximum driving force observed between calls to PrintDiagnostics()
        bool Averaging;                                                             ///< 控制参数，用于平均界面上的驱动力。Control parameter for averaging the driving force over the interface.
        std::string EnDensUnits;                                                    ///< VTK输出的能量密度单位Energy density Units for VTK output
        double CutOff;                                                              ///< 驱动力截断值Driving force cutoff value

        private:
        double maxPsi;
    };
}//namespace opensim

#endif