

#ifndef CHEMICALPROPERTIES_H
#define CHEMICALPROPERTIES_H

#include "Tools/Includes.h"
#include "Chemistry/ThermodynamicPhase.h"

namespace opensim
{

// Storages needed for TracerDiffusion
struct delta                                                                  ///  Structure of in- and outgoing composition
{
	/**This structure stores all incoming and outgoing tracer fluxes, used in the
	diffusion scheme*/
	double in;
	double out;
	delta& operator*= (double value)
	{
		in *= value;
		out *= value;
		return *this;
	}
	delta& operator+= (delta& locDelta)
	{
		in += locDelta.in;
		out += locDelta.out;
		return *this;
	}
	delta& operator-= (delta& locDelta)
	{
		in -= locDelta.in;
		out -= locDelta.out;
		return *this;
	}
};

class PhaseField;

/***/
class ChemicalProperties : public OPObject                                      ///  储存体系的热力学性质类的定义
{
 public:
    //Initialize + ReadInput
    ChemicalProperties(std::vector<std::string> ELnames,
                       std::vector<std::string> PHnames)
     {
        Ncomp = int(ELnames.size());
        Nphases = int(PHnames.size());
     }
    ChemicalProperties(){}
    void Initialize(Settings& locSettings, const int boundarysize);             ///< 通过读取之前输入维度来初始化储存空间
    void InitializeTracer(const int boundarysize);								/// 初始化追踪器
	void ReadInput(Settings& locSettings, std::string FileName);                                       ///< 从输入文件中读取输入数据
//   void ReadManualInput(std::string FileName);                                 ///< -废弃的-
    void ReadInputTracer(std::fstream& inp , int moduleLocation);				/// 从追踪器中读取输入参数
    void ReadCompositionInput(std::string FileName);                            ///< 从每个相中读取初始成分
    void Read(BoundaryConditions& BC, PhaseField& Phi,
              int tStep, bool legacy_format = true);                            ///< 从文件中读取原始成分
    //Set Boundary Conditions and Composition limits
    void SetBoundaryConditions(BoundaryConditions& BC);                         ///< 设置边界条件
    void SetBoundaryConditionsTracer(BoundaryConditions& BC);                         ///< 设置边界条件
    void SetLimitsBoundaryConditions(BoundaryConditions& BC);                   ///< 为极限设置边界条件 //TODO: obsolete
    //Large Deformation + Moving Frame Functions
    void MoveFrame(int di, int dj, int dk, BoundaryConditions& BC);             ///< 对应X，Y，Z将数据转换到 di, dj and dk的储存方式 
    void Remesh(int newNx, int newNy, int newNz, BoundaryConditions& BC);       ///< 在模拟盒子改变维度的时候，调整所有储存的维度
    //Write Output
    void Write(int tStep, bool legacy_format = true);                           ///< 将原始成分写进一个文件
    void WriteVTK(PhaseField& Phase, int tStep);                                ///< 将成分写为VTK文件格式 (.vtk file)
    void WriteStatistics(Composition& Cx, PhaseField& Phi, int tStep, double dt);///< 将成分统计写入文件, 输入：时间步 //TODO: obsolete, will be deleted!
    //Calculate Mole/Weight Fractions and total Moles from Site Fractions
    void CalculateMoleFractions(Composition& Cx);                               ///< 计算所有相的摩尔分数
    double Nmoles(Composition& Cx, PhaseField& Phi, int x, int y, int z);       ///< 返回计算后每个点的摩尔数
    //3D Storages for Phases, Components and compositions
    std::vector<double> getTotalComposition(Composition& Cx, PhaseField& Phi,int x,int y,int z); ///< 计算所有组分的质量百分比
    std::vector<double> getWeightPercent(Composition& Cx, PhaseField& Phi, int x, int y, int z); ///< 计算所有组分的质量百分比
    void SetTracerComposition();
    void GetInitialTracerAmount();
    double GetTracerAmount(int trac);
    double GetTracerFraction(Composition &Cx, int x, int y, int z, int trac);
    std::vector<Element> Component;                                             ///< 列出目前相中所有的化学组分
    std::vector<ThermodynamicPhase> Phase;                                      ///< 在模拟过程中考虑的相
    //Static Storages, set at Initialization
    int Nx, Ny, Nz;                                                             ///< 在x，y，z方向计算区域的大小
    double Threshold;                                                           ///< 保持扩散流的最小相分数
    std::vector<double> TotInitial;                                             ///< 每个组分的初始量
    std::vector<std::string> Phasenames;                                        ///< 子集挑选出的相名称
    std::vector<std::string> Elementnames;                                      ///< 子集挑选出的相名称
    std::vector<std::string> RefElements;
    //函数，访问存储空间
    int PhaseNumber(std::string name);                                          ///< 返回热力学相的索引数
    //TODO:以下函数的私有化:
    int Ncomp;                                                                  ///< 返回组分的数目
    int Nphases;                                                                ///< 返回相的数目
    int TotalNsubs;                                                             ///< 返回一个整数来表示系统中的亚点阵数
    std::vector<int> Ncons;                                                     ///< 返回一个向量来表示每个相的组成数
    std::vector<int> Nsubs;                                                     ///< 返回每个阶段的子格子数量的向量
    std::vector<double> Nsites;                                                 ///< 返回每个阶段的站点总数的矢量
    std::vector<std::vector<int> > SublIdx;                                     ///< 返回一个矩阵，其中包含每个子格的开始和结束成分
    std::vector<std::vector<int> > ConsIdx;                                     ///< 返回一个包含每个组成部分的索引的矩阵
    std::vector<std::vector<double> > NsitesWithI;                              ///< 返回一个矩阵
    //Temporary Storage, might be moved elsewhere
    std::vector<std::string> PhasenamesTQ;                                      ///< TQ接口中阶段的名称
    std::vector<std::string> ElementnamesTQ;                                    ///< TQ接口中元素的名称

    double Eta;                                                                 ///< 界面宽度
    double dx, dt;
    double TotalMolarVolume;
    void CalculateTotalMolarVolume(Composition& Cx, PhaseField& Phi);
    void CalculateTotalAverage(Composition& Cx, PhaseField& Phi);
    std::vector<double> TotalAverage;
    std::vector<int> CalculateSortedElementMatrix(void);
    std::vector<int> SortedElementMatrix;
    bool ElementsAreSorted;


    bool Isotope; 																///< 处理模拟中的同位素

    std::vector<int> IsotopeDefined;                                            ///< 包含有示踪原子的那些成分的索引
    int NumIsotopes;
    std::vector<double> InitialTracerAmount;                                    ///< 包含有同位素的所有元素
    std::vector<double> InitialTracerPosition;                                    ///< 包含同位素的所有位置
    Storage3D<double,1> TotalTracer;                                            ///< 总的示踪物组成存储
    Storage3D<delta,4> PhaseDotTracer;                                          ///< 总的示踪物组成存储
    std::vector<double> minitialTracer;

    enum class DiffusionModel {None, EquilibriumPartitioning, FiniteInterfaceDissipation };
    DiffusionModel diffMod;                                                     ///< 扩散模型的类型

 protected:
    int AtStart;                                                                ///< 表明第一个时间步

 private:
    int getNcomp(void);                                                         ///< 返回组元的数量
    int getNphases(void);                                                       ///< 返回阶段的数目
    int getTotalNsubs(void);                                                    ///< 返回系统中子格总数的整数
    std::vector<int> getNcons(void);                                            ///< 返回每个相的组分数量的矢量
    std::vector<int> getNsubs(void);                                            ///< 返回每个阶段的子格子数量的向量
    std::vector<double> getNsites(void);                                        ///< 返回每个阶段的站点总数的矢量
    std::vector<std::vector<int> > getSublIdx(void);                            ///< 返回一个矩阵，其中包含每个子格的开始和结束成分
    std::vector<std::vector<int> > getConsIdx(void);                            ///< 返回一个包含每个组成部分的索引的矩阵
    std::vector<std::vector<double> > getNsitesWithI(void);                     ///< 返回一个矩阵
    void PrintData(void);                                                       ///< 打印化学数据到屏幕
    void CalculateMolefractionLimits(void);                                     ///< 计算每个相的摩尔分数限制
    double GetMolarVolume(std::string Element);                                 ///< 返回给定元素的摩尔体积
    double GetAtomicWeight(std::string Element);                                ///< 返回给定元素的原子量
};
}
#endif//CHEMICALPROPERTIES_H

