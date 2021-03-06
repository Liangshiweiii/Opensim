

#ifndef SYSTEMDIMENSIONS_H
#define SYSTEMDIMENSIONS_H

#include "Tools/Includes.h"

namespace opensim
{
//class OPObject;

class SystemDimensions /*: public OPObject*/                                        ///< 存储体系爱三维尺度以及它们随时间的改变，如果体系形变之后自动重新分配网格。 
{
 public:
    SystemDimensions(void)                                                      ///< Default constructor默认构造函数
    {
        //Initialize();
        //ReadInput();
    }
    
    SystemDimensions(std::string InputFileName)                                 ///< 使用用户特定的输入文件来初始化
    {
        Initialize();
        ReadInput(InputFileName);
    }
     
    std::string thisclassname;                                                  ///< 对象的名称
    std::string LengthUnits;                                                    ///< 长度单位
    
    int Nx;                                                                     ///< X方向的实际格点数
    int Ny;                                                                     ///< Y方向的实际格点数
    int Nz;                                                                     ///< Z方向的实际格点数

    int dNx;																	///< X方向为了减小维度模拟的模板步长
    int dNy;																	///< Y方向为了减小维度模拟的模板步长
    int dNz;																	///< Z方向为了减小维度模拟的模板步长

    int newNx;                                                                  ///< X方向上新的格点数
    int newNy;                                                                  ///< Y方向上新的格点数
    int newNz;                                                                  ///< Z方向上新的格点数 
    
    int maxNx;                                                                  ///< X方向上最大的格点数 
    int maxNy;                                                                  ///< Y方向上最大的格点数 
    int maxNz;                                                                  ///< Z方向上最大的格点数 
    
    int OffsetX;                                                                ///< X方向的调节量
    int OffsetY;                                                                ///< Y方向的调节量
    int OffsetZ;                                                                ///< Z方向的调节量
    
    int TotalX;                                                                 ///< X方向的总的系统大小 (在MPI parallelism例子中)
    int TotalY;                                                                 ///< Y方向的总的系统大小 (在MPI parallelism例子中)
    int TotalZ;                                                                 ///< Z方向的总的系统大小 (在MPI parallelism例子中)
    
    double dx;                                                                  ///< 空间步长
    double iWidth;																///< 用网格大小表示的界面宽度
    
    bool RemeshingAllowed;                                                      ///< 允许重排网格？ (Yes or No)  
    bool initialized = false;
    void Initialize(void);                                                      ///< 初始化此类
    void ReadInput();                                                           ///< 从默认输入文件中读取输入
    void ReadInput(std::string InputFileName);                                  ///< 从特定输入文件中读取输入
    void Resize(int size_x, int size_y, int size_z, int tStep);                 ///< 在新的维度重新标定系统大小
    void Read(int tStep);                                                       ///< 从文件中读取历史系统维度
    void Write(int tStep);                                                      ///< 将历史系统维度写入文件 
    void AddForRemesh(OPObject& obj);                                           ///< 往ObjectsToRemesh类中加入对象
    void RemeshAll(void);                                                       ///< 调用Remesh()作用于ObjectsToRemesh的所有对象
    void Remesh(std::string ObjectName, int size_x, int size_y, int size_z);    ///< 调用Remesh()作用于object(s)以特定命名
    std::vector<std::array<int, 4>> History;                                    ///< 根据模拟时间来存储维度演化历史，语法: {{tStep1, Nx1, Ny1, Nz1}, ... ,{tStepN, NxN, NyN, NzN}}
    std::vector<OPObject*> ObjectsToRemesh;                                     ///< Stores pointers to objects sensitive to remeshing
    std::vector<std::string> NamesOfObjectsToRemesh;                            ///< Stores the names of objects sensitive to remeshing
    SystemDimensions& operator= (const SystemDimensions& rhs);                  ///< 复制运算符
};

}// namespace openphase

#endif
