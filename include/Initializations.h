#ifndef INITIALIZATIONS_H
#define INITTALIzATIONS_H

#include "Tools/Includes.h"

namespace opensim{
    class BoundaryConditions;
    class ElasticProperties;
    class Orientations;
    class PhaseField;
    class Settings;

    class Initializations:public OPObject{
        public:
        static void CreateMicrostructure(std::string InputFileName, PhaseField& Phase, BoundaryConditions& BC, Settings& OPSettings);//微结构创建函数声明

	    static std::vector<std::vector<int>> QuasiRandomNuclei(PhaseField& Phase, Settings& locSettings, int phaseIndex, int dist, int seed = -1, int offset = 2);//随机种子函数的一种
	    static void RandomNuclei(PhaseField& Phase, Settings& locSettings, int phaseIndex, int Nparticles, int seed = -1);//随机种子函数的声明
        static int Ellipsoid(PhaseField& Phase, int PhaseIndex, double RadiusX,
            double RadiusY, double RadiusZ, double x0, double y0, double z0,
            BoundaryConditions& BC, Settings& locSettings);           //椭圆相区域的设置函数声明
        static int SectionalPlane(PhaseField& Phase, const int PhaseIndex,
            const dVector3 Point, const dVector3 Orientation,
            const BoundaryConditions& BC, const Settings& locSettings,
            const bool NewGrain = true, const int FieldIndex = 0);             /// 在已有平面背面的初始化另一个平面
        static int Single(PhaseField& Phase, int PhaseIndex, BoundaryConditions& BC,
                Settings& locSettings);                                               ///单相的初始化函数的声明
        static int Zlayer(PhaseField& Phase, int PhaseIndex, int Position, 
                        int Thickness, BoundaryConditions& BC, Settings& locSettings);    ///Z平面设置函数的声明
        static int Sphere(PhaseField& Phase, int PhaseIndex, double Radius,
                double x0, double y0, double z0, BoundaryConditions& BC,
                Settings& locSettings);                                             ///球型相区域设置函数的声明
        static std::vector<int> Fractional(PhaseField& Phase,
                int MajorityPhaseIndex, int MinorityPhaseIndex,
                double MinorityPhaseLayerThickness, BoundaryConditions& BC,
                Settings& locSettings);                                           ///单相分配函数的声明
        static std::vector<int> FractionalSharpInterface(PhaseField& Phase,
                int MajorityPhaseIndex, int MinorityPhaseIndex,
                double MinorityPhaseLayerThickness, BoundaryConditions& BC,
                Settings& locSettings);                                           ///尖锐界面单相分配函数的声明
        static std::vector<int> ThreeFractionals(PhaseField& Phase,
                int MajorityPhaseIndex, double MajorityPhaseLayerThickness,
                int MinorityPhaseIndex1, double MinorityPhaseLayerThickness1,
                int MinorityPhaseIndex2, BoundaryConditions& BC,
                Settings& locSettings);                                           ///三相分配函数的定义
        static std::vector<int> TwoDifferentWalls(PhaseField& Phase,
                int ChannelPhaseIndex, int WallsPhaseIndex, double WallsThickness,
                BoundaryConditions& BC, Settings& locSettings);
        static std::vector<int> TwoWalls(PhaseField& Phase, int ChannelPhaseIndex,
                int WallsPhaseIndex, double WallsThickness, BoundaryConditions& BC,
                Settings& locSettings);
        static void Cube(PhaseField& Phase, int PFIndex, double lx, double ly,
                double lz, double x0, double y0, double z0, BoundaryConditions& BC,
                Settings& locSettings);                                           ///立方结构设置函数的声明
        static int CubeSimple(PhaseField& Phase, const int PhaseIndex,
                const double Lx, const double Ly, const double Lz,
                const int i0, const int j0, const int k0,
                const BoundaryConditions& BC, const Settings& locSettings);       ///立方简单结构设置函数的声明
        static void Cylinder(PhaseField& Phase, int PhaseFieldIndex, double Radius,
                double length, int Axis, double x0, double y0, double z0,
                Settings& locSettings);                                           ///圆柱相结构设置函数的声明

        // Special initialization functions:
        static int  FillGrainWithRandomSpheres(PhaseField& Phase,
                int ParrentPhaseFieldIndex, int SpheresPhaseIndex, double x0,
                double y0, double z0, double Lx, double Ly, double Lz, double MinR,
                double MaxR, BoundaryConditions& BC, Settings& locSettings);      ///晶粒中随机球状形核函数的声明
        static int  FillGrainWithSpheres(PhaseField& Phase,
                int ParrentPhaseFieldIndex, int SpheresPhaseIndex,
                double SphereRadius, BoundaryConditions& BC, Settings& locSettings);///晶粒汇总球状形核函数的声明
        static int  FillGrainWithSpheres(PhaseField& Phase,
                int ParrentPhaseFieldIndex, int SpheresPhaseIndex, double x0,
                double y0, double z0, double Lx, double Ly, double Lz,
                double SphereRadius, BoundaryConditions& BC, Settings& locSettings);///晶粒汇总球状形核函数的重载
        static int DiffusionCouple(PhaseField& Phase, int MajorityPhaseIndex,
                int MinorityPhaseIndex, double MinorityPhaseLayerThickness,
                BoundaryConditions& BC, Settings& locSettings);                     ///扩散偶函数的声明
        /*static int MechanicalWorkBench(PhaseField& Phase, int nrVoidCells,
                int voidPhIDX, int matrixPhIDX, BoundaryConditions& BC,
                Settings& locSettings);*/
        static int SphereInGrain(PhaseField& Phase, int ParrentPhaseFieldIndex,
                int PhaseIndex, double Radius, double x0, double y0, double z0,
                BoundaryConditions& BC, Settings& locSettings);
        static int TwoDimEBSD(std::string filename, std::vector<int> columns,
                PhaseField& Phase, BoundaryConditions& BC, Settings& locSettings);  /// Reads microstructure from EBSD result files从EBSD结果中读取微结构
        static int TwoDimEBSDWithOrientations(std::string filename,
                std::vector<int> columns, std::string anglerepresentation,
                PhaseField& Phase, ElasticProperties& EP, Orientations& OR,
                BoundaryConditions& BC, Settings& locSettings);                     /// Reads microstructure from EBSD result files从EBSD结果中读取微结构211
        static void VoronoiTesselation(PhaseField& Phase, BoundaryConditions& BC,
                Settings& OPSettings, int Ngrains, int MatrixPhase);                                         /// 使用Voro++库来生成Voronoi晶粒结构图
        static void Young3(PhaseField& Phase, int PhaseIndex,
                BoundaryConditions& BC, Settings& locSettings);
        static void Young4(PhaseField& Phase, int PhaseIndex,
                BoundaryConditions& BC, Settings& locSettings);
        static void Young4Periodic(PhaseField& Phase, int PhaseIndex,
                BoundaryConditions& BC, Settings& locSettings);

        static void Read(PhaseField& Phase, std::string FileInputName,
                BoundaryConditions& BC, Settings& locSettings);                                   ///读取相信息函数的声明  
        static void ReadSubset(PhaseField& Phase,    
                std::string PFFileName, 
                std::string GSFileName,
                int offsetX,
                int offsetY,
                int offsetZ,
                const BoundaryConditions& BC,
                int inpNx = 1,
                int inpNy = 1,
                int inpNz = 1);                                                     ///< 从可能不同(更大)大小的文件中读取原始(二进制)相位字段的子集
        static void ReadSubsetGeneral(PhaseField& Phase, 
                std::string PFFileName,                                             ///< 相位字段输入文件名
                std::string GSFileName,                                             ///< 晶界统计输入文件名
                int offsetInpX,                                                     ///< 输入读数的起始点
                int offsetInpY,
                int offsetInpZ,
                int sizeInpX,                                                       ///< 从offsetInpX/Y/Z开始读取输入数据的大小
                int sizeInpY,
                int sizeInpZ,
                int totalSizeInpNx,                                                 ///< 关于遗留格式，输入数据的总大小
                int totalSizeInpNy,
                int totalSizeInpNz,
                int offsetLocalX,                                                   ///< 粘贴到本地数据的起点
                int offsetLocalY,
                int offsetLocalZ,
                std::initializer_list<bool> newGrain,                               ///< 创建新的粒度为true，复制到目标中的相同粒度为false，如果没有指定值，则假设为true
                const BoundaryConditions& BC);                                      ///< 从可能不同大小的文件中读取原始(二进制)相位字段的子集，粘贴到自定义位置

        private:
    };
}

#endif


