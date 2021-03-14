#ifndef SETTINGS_H
#define SETTINGS_H

#include "Tools/Includes.h"
#include "SystemDimensions.h"
#include "Chemistry/ChemicalProperties.h"



namespace opensim
{
/**Directory Structure Settings:*/

class Settings                                                                  ///< OpenPhase 设置模块. 读取储存模拟的所有变量 
{
 public:
    std::string thisclassname;                                                  ///< Object's class name定义模块设置类名
    //std::string DefaultInputFileName;                                         ///< Default input file name默认输入的文件名
    std::string UnitsOfLength;                                                  ///< Units of length长度单位
    std::string UnitsOfMass;                                                    ///< Units of mass质量单位
    std::string UnitsOfTime;                                                    ///< Units of time时间单位
    std::string UnitsOfEnergy;                                                  ///< Units of energy能量单位

    int Nx;                                                                     ///< Number of grid points in X directionX在X方向的格点数
    int Ny;                                                                     ///< Number of grid points in Y direction在Y方向的格点数
    int Nz;                                                                     ///< Number of grid points in Z direction在Z方向的格点数

// New
    SystemDimensions SD;                                                ///< System dimensions object创建系统维度类和化学性质类
	ChemicalProperties Cp;
// end New

    int Nphases;                                                                ///< Number of phases per grain每个晶粒中的相数
    int Ncomp;                                                                  ///< Number of the different components不同的组元数
    double dt;                                                                  ///< Initial time step时间步长
    int nSteps;                                                                 ///< Number of time steps时间步长的个数
    int tFileWrite;                                                             ///< Writing distance to the file (in time steps)输出到文件的时间步
    int tScreenWrite;                                                           ///< Writing distance to the screen (in time steps)输出到屏幕的时间步
    int tStart;                                                                 ///< Starting time value (in time steps). Used for settig the restart time step起始时间步
    int Restart;                                                                ///< Restart switch ("No" if restart is OFF, "Yes" if ON)重新开始阈值
    int tRestartWrite;                                                          ///< Restart distance (in time steps)重启时间步

    int nOMP;                                                                   ///< Number of OpenMP Threads线程数


    double dx;                                                                  ///< Grid spacing空间步长
    double iWidth;                                                              ///< Interface width in grid points界面宽度（格点单位）
    double Eta;                                                                 ///< Interface width in units of length界面宽度（长度单位）

    std::string SimulationTitle;

    void Initialize(void);
    void ReadInput();
    void ReadInput(std::string InputFileName);

    Settings& operator= (const Settings& rhs);                                  ///< Copy operator for Settings class设置类的复制运算符
};

}// namespace opensim

#endif