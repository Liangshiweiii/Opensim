

#ifndef PHASEFIELD_H
#define PHASEFIELD_H

#include "Tools/Includes.h"
#include "Tools/UserInterface.h"
#include "Settings.h"
#include "Tools/Node.h"
#include "Tools/NodeV.h"
#include "Tools/Tensor.h"
#include "GrainInfo.h"

namespace opensim
{

class ElasticitySolverSpectral;
class BoundaryConditions;
class ChemicalProperties;
class Velocities;
class FlowSolverLBM;

/******************************************************************************/
class PhaseField : public OPObject                                              ///< 相场类，它存储相场信息，并对它们执行基本操作。
{
 public:

    PhaseField(){ };
    ~PhaseField()
    {
#ifdef _OPENMP
        destroy_omp_locks();
#endif
    };
    PhaseField(const Settings& locSettings, unsigned int boundary = 2);

    int Initialize(const Settings& locSettings, unsigned int boundary = 2, int df = 1);    ///< 初始化驱动力类的存储和初始变量。
    int  PlantGrainNucleus(int PhaseIndex, int x, int y, int z);                ///< 在“PhaseIndex”相的(x,y,z)位置放置晶核，返回得到的相场指数
    int  AddGrainInfo(int PhaseIndex);                                          ///< 为“PhaseIndex”相位添加新的粒度信息，返回得到的相位字段索引
    void Clear();                                                               ///< 清空相场储存空间
    void CalculateLaplacians(void);                                             ///< 计算拉普拉斯算子并将其存储在相场节点中
    //Tensor<double,1> Fractions(const int i, const int j, const int k) const;  ///< 返回给定位置的相场分数
    NodeV Gradients(const int i, const int j, const int k) const;               ///< 返回给定位置上所有相场的梯度
    NodeV Normals(const int i, const int j, const int k) const;                 ///< 返回给定位置的所有相场对的界面法线
    NodeV NormalsPhase(const int i, const int j, const int k) const;            ///< 返回给定位置的所有相场对的界面法线
    std::pair<Node,Node> PrincipleCurvatures(const int i, const int j,
        const int k) const;                                                     ///< 返回给定位置所有相场的主曲率

    void CalculateVolumes();                                                    ///< 收集每个相场的体积。
    void CalculateVolumesStepOne();                                             ///< 收集每个相场的体积。
    void CalculateVolumesStepTwo();                                             ///< 收集每个相场的体积。

    void MergeIncrements(const BoundaryConditions& BC, const double dt, const bool finalize = true);///< 将增量合并到相场中
    void NormalizeIncrements(const BoundaryConditions& BC, const double dt);    ///< 规范化界面场，以使合并后得到的相场不会逸出间隔[0,1]且总和为1。

    void Read(const BoundaryConditions& BC, int tStep);                         ///< 从文件中读取原始（二进制）相场
    void Read(std::string FileName, const BoundaryConditions& BC);              ///< 从文件中读取原始（二进制）相场
    void Finalize(const BoundaryConditions& BC, const bool finalize = true);    ///< 在此时间步中结束相场计算
    void MoveFrame(const int dx, const int dy, const int dz,
                   const BoundaryConditions& BC);                               ///< 分别在x，y和z方向上将存储中的数据移动dx，dy和dz（它们应为0，-1或+1）。
    void ConsumePlane(const int dx, const int dy, const int dz,
                   const int x, const int y, const int z,
                   const BoundaryConditions& BC);                               ///分别在x，y和z方向上将存储中的数据移动dx，dy和dz（它们应为0，-1或+1）。
    void Advect(Velocities& Vel, BoundaryConditions& BC, double dt, int scheme = Upwind);///< 提出相场
	void Advect(Velocities& Vel, BoundaryConditions& BC, FlowSolverLBM& LBM, double dt, int scheme = Upwind);///< 提出相场
    void PrintPFVolumes() const;                                                ///< 将所有相场的体积打印出来
    void PrintPointStatistics(const int x, const int y, const int z) const;     ///< 将本地显示的相场及其值打印到屏幕上。
    void PrintVolumeFractions(std::vector<std::string> Names);                  ///< 打印所有热力学相的体积分数
    void PrintVolumeFractions(ChemicalProperties& CP);                          ///< 打印所有热力学相的体积分数
    void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC); ///< 在保持数据不变的情况下改变网格
    void SetBoundaryConditions(const BoundaryConditions& BC);                   ///< 设置边界条件
    void SetIncrementsBoundaryConditions(const BoundaryConditions& BC);         ///< 为相场增量设置边界条件
    void SetFlags();                                                            ///< 设置标记界面的标志
    void Write(int tStep);                                                      ///< 将原始相场信息写入文件
    void WriteAverageVolume(const int tStep, const int PhaseIndex) const;       ///< 将给定相场的平均体积信息写入文件
    void WriteDistortedVTK(int tStep, ElasticitySolverSpectral& ES, double scale);     ///< 将扭曲网格上的相场写入到VTK文件中
    void WriteGrainsStatistics(const int tStep);                                ///< 在现在的时间步将晶粒和交界处的数据写入文件中
    void WriteIndividualPhaseFieldValuesVTK(const int tStep,
            const std::initializer_list<int> FieldIndices) const;               ///< 将单独的相场值写入VTK文件中
    void WriteLaplaciansVTK(int tStep, int PhiIndex);                           ///< 将给定相场的拉普拉斯算子以VTK格式写入文件
    void WriteVTK(const int tStep, const bool CurvatureOutput = false) const;   ///< 以VTK的格式将相场写入文件
    void WriteVTKData(std::stringstream& buffer,
            const bool CurvatureOutput = false) const;                          ///< 将相场以VTK格式写入缓冲区，以便可以将多个VTK输出放入一个文件中

//    void CalculateGrainboundaryPotentialContrib(ThermodynamicFunctions& TF,
//                                                     Settings& locSettings);  ///< 计算晶界对扩散势的贡献
    bool PhaseFieldPresent(const int i, const int j, const int k,
            const int Index) const;                                             ///< 如果（i，j，k）存在相场，则返回true

    std::vector<int> GetPresentPhaseFields() const;                             ///< 返回带有当前相场索引的有序列表
    std::vector<int> ReturnVicinityPhaseFields(const int i,
            const int j, const int k) const;                                    ///< 返回点（i，j，k）和相邻节点处的相场列表

    bool Interface(const int x, const int y, const int z) const                 ///< 指示接口的位置。
    {
        return (Fields(x,y,z).flag > 1);
    };
    Storage3D<   Node, 0 >  Fields;                                             ///< 相场存储
    Storage3D<   Node, 0 >  FieldsDot;                                          ///< 相场增量存储
	Storage3D<double, 1> Fractions;                                             ///< 每个点的热力学相分数

	void CalculateFractions();				///< 返回热力学相分数

    double dx;                                                                  ///< 网格间距
    int    Nx;                                                                  ///< 系统的X维度
    int    Ny;                                                                  ///< 系统的Y维度
    int    Nz;                                                                  ///< 系统的Z维度
    int    Nphases;                                                             ///< 热力学相数
    double Eta;                                                                 ///< 物理单位的接口宽度
    double iWidth;                                                              ///< 界面宽度（以网格点为单位）

    bool NucleationPresent;                                                     ///< 如果存在任何相的核，则为true，否则为false
    double RefVolume;                                                           ///< 参考体积
    GrainInfo FieldsStatistics;                                                 ///< 相场统计。包含有关每个相场的位置，速度，方向等的信息

    PhaseField& operator= (const PhaseField& rhs);                              ///< PhaseField类的复制运算符

    void CombinePhaseFields(const int PhaseIndex);                              ///< 使用索引PhaseIndex将同一相的相场合并到相场中
    void SelectiveCombinePhaseFields(BoundaryConditions& BC,                    ///<将sourceIndex的相场合并到索引为targetPhaseIndex的相场
            const int SourcePhaseIndex, const int TargetPhaseIndex);
    std::vector<int> GetMaxPhaseFieldOverlap(const int thPhase1,                 ///< Returns a vector containing the phaseindices thPhase1 and thPhase2 that have the biggest overlap, also the number of interface points between the phasefields. vector contains 1. phasefield of thPhase1, 2. phasefield of thePhase2, 3. number of overlap points. All values are -1 if no overlap is found.
            const int thPhase2);
    Matrix<int> GetPhaseFieldOverlap(const int thPhase1,                         ///< 返回一个矩阵，其中包含在矩阵的列和行中表示的相索引的重叠。仅设置对角线以上的值。
            const int thPhase2);
    int resolution;                                                             ///< 相场的网格分辨率。可以是单或双
    void set_omp_lock(int i, int j, int k)
    {
#ifdef _OPENMP
        omp_set_lock(&omp_locks(i, j, k));
#endif //_OPENMP
    }
    void unset_omp_lock(int i, int j, int k)
    {
#ifdef _OPENMP
        omp_unset_lock(&omp_locks(i, j, k));
#endif //_OPENMP
    }
    void init_omp_locks();
    void destroy_omp_locks();

 protected:
 private:
#ifdef _OPENMP
    Storage3D<omp_lock_t, 0> omp_locks;                                             ///< omp lock storage
#endif
};

}// namespace opensim
#endif

