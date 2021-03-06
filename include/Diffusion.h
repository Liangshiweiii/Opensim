#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "Tools/Includes.h"
#include "Diffusion/DiffusionStencil.h"

namespace opensim{
    class Settings;
    class PhaseField;
    class DrivingForce;
    class Composition;
    class Temperature;
    class InterfaceMobility;
    class BoundaryConditions;

    class Diffusion : public OPObject{
        public:
        Diffusion(){}
        Diffusion(const Settings& locSettings,
                const int boundary = 1);                                            ///<  初始全局设定Initializes global settings
        void Initialize(const Settings& locSettings,
                const int boundary = 1);                                            ///<  初始全局设定Initializes global settings
        using OPObject::ReadInput;
        void ReadInput(const std::string InputFileName);                            ///<  从一个文件中读取输入参数Reads input parameters from a file
        void GetDrivingForce(PhaseField& Phase, Composition& Cx, Temperature& Tx,
                                                                DrivingForce& dGab);///<  计算每个点的驱动力Calculates the driving force for each point
        void GetDrivingForceLoop(PhaseField& Phase, Composition& Cx, Temperature& Tx,
                                                                DrivingForce& dGab);///<  计算每个点的驱动力Calculates the driving force for each point
        double GetDrivingForceAlpha(PhaseField& Phase, Composition& Elements,
                                                    Temperature& Tx, int alpha,
                                                            int i, int j, int k);  ///<  计算已有点的alpha相的驱动力Calculates the driving force for phase alpha in a given point
        double SetInitialComposition(PhaseField& Phase, Composition& Cx);           ///<  设置初始成分Sets the initial composition

        void CalculateLocalPhaseConcentrations(PhaseField& Phase, Temperature& Tx,
                                            Composition& Cx, int i, int j, int k);///<  将总成分分配到给定点的每个相的成分中
        void CalculatePhaseConcentrations(PhaseField& Phase,
                                                Temperature& Tx, Composition& Cx);///<  将总成分分配到给定点的每个相的成分中
        void SetDiffusionCoefficients(PhaseField& Phase, Temperature& Tx);          ///<  设置每个点的扩散系数
        void Solve(PhaseField& Phase, Composition& Cx,
                                Temperature& Tx, BoundaryConditions& BC, double dt,
                                bool InternalDiffusivities = true);                 ///<  考虑到交叉项，可以一步一步计算总浓度的变化
        void CalculateInterfaceMobility(PhaseField& Phase, Composition& Cx,
                            BoundaryConditions& BC, InterfaceMobility& IMobility);///<  计算依赖于浓度的迁移率
        void RestoreStoichiometric(PhaseField& Phase, Temperature& Tx,
                                                                Composition& Cx);///<  恢复完全生长的化学计量池中的精确化学计量组成
        void RestoreStoichiometricThreadSafe(PhaseField& Phase, Temperature& Tx,
                                                                Composition& Cx);///<  恢复完全生长的化学计量池中的精确化学计量组成
        void Remesh(int newNx, int newNy, int newNz);                               ///<  在保留数据的同时重新划分存储空间
        void CalculateDiffusionIncrements(PhaseField& Phase, Composition& Elements,
                                                                Temperature& Tx);///<  计算菲克的扩散浓度增量
        void CalculateAntitrappingIncrements(PhaseField& Phase,
                                                            Composition& Elements);///<  Calculates antitrapping concetration increments
        void LimitAntitrappingIncrements(PhaseField& Phase,
                                            Composition& Elements, Temperature& Tx);///<  Limits antitrapping concetration increments
        void LimitAntitrappingIncrementsThreadSafe(PhaseField& Phase,
                                            Composition& Elements, Temperature& Tx);///<  Limits antitrapping concetration increments
        void CalculateLimits(PhaseField& Phase, Composition& Elements, double dt);  ///<  计算的增量极限，以使浓度保持在物理范围内
        void LimitDiffusionIncrements(PhaseField& Phase,
                                            Composition& Elements, Temperature& Tx);///<  限制Fick的扩散速度增量
        void ApplyIncrements(PhaseField& Phase, Composition& Elements, double dt);  ///<  应用限制的增量
        double ReportMaximumTimestep(PhaseField& Phase);
        void SetPhaseFractions(PhaseField& Phase);                                  ///<  设置每个点的热力学相分数
        void PrintPointStatistics(int x, int y, int z);                             ///<  打印点中的扩散系数
        int Nx;                                                                     ///<  沿X的内部计算域的大小
        int Ny;                                                                     ///<  沿Y的内部计算域的大小
        int Nz;                                                                     ///<  沿Z的内部计算域的大小

        bool LimitingNeeded;
        int Comp;                                                                   ///<  给定实例的扩散求解器的化学成分索引
        int Nphases;                                                                ///<  热力学相数
        double Precision;                                                           ///<  成分评估的精度
        double Eta;                                                                 ///<  界面宽度
        double dx;                                                                  ///<  网格宽度
        double dx_1;                                                                ///<  网格间距的倒数
        double dx_2;                                                                ///<  网格间距的平方的倒数
        double R;                                                                   ///<  通用气体常数
        double TotalMass;                                                           ///<  给定求解器实例处理的组件总数

        std::vector<int> Stoichiometric;                                            ///<  每个热力学阶段的化学计量标志
        std::vector<double> DC0;                                                    ///<  每个热力学相的溶质扩散系数
        std::vector<double> AE;                                                     ///<  每个热力学阶段的扩散激活能
        Tensor<double, 2> Pt;                                                       ///<  线性相图的成对分配系数
        Tensor<double, 2> mL;                                                       ///<  线性相图的成对液相线斜率
        Tensor<double, 2> mL_1;                                                     ///<  线性相图的成对液相线斜率的倒数
        Tensor<double, 2> Ts;                                                       ///<  线性相图的液相线-固相线或固溶线相交的成对温度
        Tensor<double, 2> Cs;                                                       ///<  液相线-固相线或固溶线相交的成对浓度线性相图

        std::vector<double> S;                                                      ///<  不同阶段的熵。它们的差异用于驱动力

        int ThereAreStoichiometricPhases;                                           ///<  指示是否存在化学计量阶段
        Storage3D<double, 1> dMu;                                                   ///<  外部对化学势的贡献。
        Storage3D<double, 1> DC;                                                    ///<  每个点的扩散系数
        Storage3D<double, 1> PhaseFractions;                                        ///<  每个点的热力学相分数

        DiffusionStencil DStencil;                                                  ///<  扩散模具。以拉普拉斯模具为基础

    };
}

#endif