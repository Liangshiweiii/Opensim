#ifndef TEMPERATURES_H
#define TEMPERATURES_H

#include "Tools/Includes.h"

namespace opensim
{
    class Settings;
    class GrainInfo;
    class PhaseField;
    class BoundaryConditions;
    class InterfaceField;
    class Velocities;

    class Temperature : public OPObject                                             ///< 温度的存储Storage for the temperature
    {
    public:
        Temperature(){};
        Temperature(const Settings& locSettings, const std::string InputFileName = "NONE");
        void Initialize(const Settings& locSettings, unsigned int boundary = 2);                             ///< 分配内存，初始化设定Allocates memory, initializes the settings
        using OPObject::ReadInput;
        void ReadInput(const std::string FileName);                                 ///< 读取初始化内部变量值Reads initial values for internal parameters.
        void MoveFrame(const int dx, const int dy, const int dz,
                    const BoundaryConditions& BC);                               ///< 将数据各自根据XYZ单独的关系构建Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in x, y or z direction correspondigly.
        void ConsumePlane(const int dx, const int dy, const int dz,
                    const int x, const int y, const int z,
                    const BoundaryConditions& BC);                               ///< 将数据按照从一个点或者到一个点的规则转化为XYZ格式Shifts the data in the storage by dx, dy and dz (they should be 0, -1 or +1) in the direction to/from (x, y, z) point.
        void Remesh(const int newNx, const int newNy, const int newNz,              ///< 在保持数据不变的情况下改变体系大小Change system size while keeping the data
                                        const BoundaryConditions& BC);
        void SetBoundaryConditions(const BoundaryConditions& BC);                   ///< 通过在重影节点中设置适当的值来设置温度的边界条件Sets boundary conditions for the temperature by setting the appropriate values in ghost nodes
        void Write(const int tStep);                                                ///< 将原始的温度写进文件中Writes the raw temperature into a file
        void Read(const int tStep);                                                 ///< 从文件中读取原始温度Reads the raw temperature from a file
        void WriteVTK(const int tStep) const;                                       ///< 以VTK格式将温度写入文件中Writes temperature in the VTK format into a file
        void WriteTemperatureGradientVTK(const int tStep) const;
        void WriteMICRESS(const int tStep) const;                                   ///< 以MiCRESS的格式将温度写入文件中Writes temperature in the MICRESS format into a file

        void PrintPointStatistics(const int x, const int y, const int z) const;     ///< 将给定的（x，y，z）点打印温度到屏幕Prints temperature at a given point (x, y, z) to screen
        void Set(const BoundaryConditions& BC, const double dt);                    ///< 根据冷速来设定温度
        void Set(const BoundaryConditions& BC, const PhaseField& Phase,
                                    const InterfaceField& Psi,
                                    const double LH, const double Cp,
                                const int PhaseIndex, const double dt);             ///< 根据冷却速度和转化产生的热量来设定温度
        void Set(const BoundaryConditions& BC, const PhaseField& Phase,
                                    const double LH, const double Cp,
                                const int PhaseIndex, const double dt);             ///< 根据冷却速度和转化产生的热量来设定温度
        void Set(const BoundaryConditions& BC, const double dTx,						///< 通过给出（+）加热/（-）冷却速率和最终温度来设置温度。
                const double Temperature, const double dt);
        void SetInitial(const BoundaryConditions& BC);                              ///< 根据起始温度和温度梯度设置初始温度
        void SetInitial(const BoundaryConditions& BC, const PhaseField& Phase,
                                                        const int PhaseIndex);      ///< 根据起始温度和温度梯度设置初始温度

        double& operator()(const int x, const int y, const int z)                   ///< 双向访问运算符
        {
            return Tx(x, y, z);
        };
        double const& operator()(const int x, const int y, const int z) const       ///< 常量引用访问运算符
        {
            return Tx(x, y, z);
        };
        double Average(void);
        Temperature& operator=(const Temperature& rhs);                             ///< 模板类的复制运算符
        void Advect(const Velocities& Vel, const BoundaryConditions& BC,
                                        const double dt, int scheme);            ///< 影响温度

        double T0;                                                                  ///< 开始温度Starting temperature
        dVector3 dT_dr;                                                             ///< 初始温度梯度Initial Temperature gradient
        dVector3 r0;                                                                ///< T0温度梯度的初始位置（相对于该位置将应用dT_dr）
        double dT_dt;                                                               ///< 冷速
        int Nx;                                                                     ///< 沿X的内部计算域的大小
        int Ny;                                                                     ///< 沿Y的内部计算域的大小
        int Nz;                                                                     ///< 沿Z的内部计算域的大小
        int Nphases;
        double dx;                                                                  ///< 空间尺寸

        Storage3D<double, 0> Tx;                                                    ///< 温度数组
        Storage3D<dVector3, 0> TxDx;                                                ///< 温度梯度
        Storage3D<double, 0> qdot;                                                  ///< 热源
        Storage3D<double, 0> TxDot;                                                 ///< 温度增量数组
    protected:
        double RefVolume;

    private:
    };
}
#endif