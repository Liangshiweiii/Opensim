#ifndef COMPOSITIONS_H
#define COMPOSITIONS_H

#include "Tools/Includes.h"
#include "Settings.h"
#include "PhaseField.h"

namespace opensim
{
    class PhaseField;
    class BoundaryConditions;
    class ElasticProperties;
    class Velocities;
    class PhaseField;

    struct CompositionDelta 
    {
        double in; 
        double out;
        
        CompositionDelta& operator*=(double m)
        {
            in  *= m;
            out *= m;
            return *this;
        }
        CompositionDelta& operator+=(CompositionDelta m)
        {
            in  += m.in;
            out += m.out;
            return *this;
        }
        CompositionDelta& operator-=(CompositionDelta m)
        {
            in  -= m.in;
            out -= m.out;
            return *this;
        }
    };

    class Composition : public OPObject                                             ///< 将组合物存储为浓度或摩尔密度等。
    {
    public:
        Composition() {}
        Composition(const Settings& locSettings, const int boundarysize = 1);
    
        void Initialize(const Settings& locSettings, const int boundarysize = 1);   ///< 初始化存储，设置内部变量。
        void InitializeEQP(const Settings& locSettings, const int boundarysize = 1);   ///< 初始化存储，设置内部变量。
        void InitializeFID(ChemicalProperties& CP, const Settings& locSettings, const int boundarysize = 1); ///< 初始化存储，设置内部变量，读取化学属性
        using OPObject::ReadInput;
        void ReadInput(std::string InputFileName);                                  ///<  从文件中读取输入参数
        void ReadInputCP(ChemicalProperties& CP, std::string InputFileName);        ///<  从文件中读取输入参数
        void ReadInputEXP(std::string InputFileName);                               ///<  从文件中读取输入参数
        void ReadInputCPFID(ChemicalProperties& CP, std::string FileName);

        Tensor<double,1> SetInitialConstituentFractions(ChemicalProperties& CP,
                                                        PhaseField& Phi,int mode);
        Tensor<double,1> SetInitialComponentFractions(ChemicalProperties& CP,
                                                    PhaseField& Phi, int mode);
        Tensor<double,1> SetInitialCompositionEQP(ChemicalProperties& CP,
                                                PhaseField& Phi);
        std::vector<double> SetInitialCompositionFID(ChemicalProperties& CP,
                                                    PhaseField& Phi);
        void CalculateTotalFractions(PhaseField& Phase);                                     ///< 计算总的成分Calculates total composition 
        void WriteVTK(int tStep);                                                   ///< 以VTK格式（.vts文件）写入成分
        void WriteVTK(PhaseField& Phi,ChemicalProperties& CP,Composition& Cx, int tStep);

        void WriteDistortedVTK(int tStep, ElasticProperties& EP, double scale);     ///< 将VTK格式的合成内容写入变形网格（.vtk文件）
        void Write(int tStep, std::string File = "Composition_" , bool legacy_format = true );                           ///< 将原始成分写入文件Writes the raw composition into a file
        void WriteSF(ChemicalProperties& CP, int tStep, bool legacy_format = true );
        void ReadSF(ChemicalProperties& CP, BoundaryConditions& BC,
                                    int tStep, bool legacy_format = true);
        void Read(ChemicalProperties& CP, BoundaryConditions& BC, int tStep,
                std::string File = "Composition_", bool legacy_format = true);    ///< 从文件中读取原始成分Reads the raw composition from a file
        void ReadSubset(BoundaryConditions& BC,
                        std::string FileName, int offsetX,
                        int offsetY,
                        int offsetZ, 
                        bool legacy_format = false,
                        int inpNx = 1, 
                        int inpNy = 1, 
                        int inpNz = 1);                                             ///<  从大小可能不同（较大）的文件中读取原始（二进制）成分的子集Read subset of a raw (binary) composition from the file of possibly different (larger) size
        void SetBoundaryConditions(const BoundaryConditions& BC);                   ///< 设置周期条件
        void SetLimitsBoundaryConditions(const BoundaryConditions& BC);             ///< 为极限设置周期条件
        void MoveFrame(const int dx, const int dy, const int dz,
                    const BoundaryConditions& BC);                               ///< 分别在x，y或z方向上将存储中的数据移动dx，dy和dz（它们应为0，-1或+1）。
        void ConsumePlane(const int dx, const int dy, const int dz, 
                        const int x, const int y, const int z,
                        const BoundaryConditions& BC);                            ///< 沿（x，y，z）点的方向将存储中的数据移动dx，dy和dz（它们应为0，-1或+1）。
        void Remesh(int newNx, int newNy, int newNz, const BoundaryConditions& BC); ///< 在保留数据的同时重新划分存储空间
        void PrintPointStatistics(int x, int y, int z);                             ///< 在给定点（x，y，z）上的成分打印到屏幕
        void WriteStatistics(int tStep, double dt);                                 ///< 将成分统计信息写入文件，输入：时间步
        void WriteStatistics(ChemicalProperties& CP, int tStep, double dt);
        
        void AdvectPhase(PhaseField& Phase, Velocities& Vel, BoundaryConditions& BC, 
                        double dt, int scheme = Upwind);                           ///< 得处相组成
        void AdvectTotal(PhaseField& Phase, Velocities& Vel, BoundaryConditions& BC, 
                        double dt, int scheme = Upwind);                           ///< 得出总成分
        void WriteVTKEQP(ChemicalProperties& CP, int tStep);
        void WriteVTKFID(PhaseField& Phi, Settings& OPSettings, int tStep);
        int Nphases;                                                                ///< 热力学相的数量Number of thermodynamic phases
        int Ncomp;                                                                  ///< 化学组分的个数Number of chemical components
        int Nx;                                                                     ///< 沿X的内部计算域的大小
        int Ny;                                                                     ///< 沿Y的内部计算域的大小
        int Nz;                                                                     ///< 沿Z的内部计算域的大小
        double dx;                                                                  ///< 网格间距
        double Threshold;                                                           ///< 仍计算扩散通量的最小相分数
        std::vector<std::string> Names;                                             ///< 相应化学成分的名称（例如Au，Cu，Na，Cl等）
        void EQPFinalize(PhaseField& Phi, ChemicalProperties& CP);



        /* template specialization integer number stands for
        * the number of extra dimensions
        * order:
        * Rank = 2: phase, component
        * Rank = 1: phase or component
        */
        Tensor<double, 2> MolarVolume;
        Tensor<double, 2> Min;                                                      ///< 每个阶段的最小成分值
        Tensor<double, 2> Max;                                                      ///< 每个阶段的最大成分值

        Storage3D<double, 2> Phase;                                                 ///< 相组成存储
        Storage3D<double, 1> Total;                                                 ///< 总成分存储

        Storage3D<CompositionDelta, 2> PhaseDot;                                    ///< 相成分增量存储
        Storage3D<double, 1> TotalDot;                                              ///< 总成分增量存储 
        Storage3D<CompositionDelta, 2> Norm;                                        ///< 储存用于相组成的归一化
        Tensor<double, 2> Initial;                                                  ///< 所有相的初始组成成分
        Tensor<double, 2> Initial2;                                                 ///< 所有相的初始组成成分
        Tensor<double, 1> Reference;                                                ///< 参考成分
        Storage3D<int, 0> Limiting;                                                 ///< 限制成分增量的存储

        std::vector<double> TotInitial;                                             ///< 每种成分的初始量

        std::vector< Storage3D< double, 1 > > ConstituentFractions;                 ///<

        Tensor<double, 1> Site2MoleFrac(ThermodynamicPhase& Phase, int i, int j, int k);
        std::vector<double> getTotalComposition(PhaseField& Phi, int x, int y, int z);
        std::vector<double> getWeightPercent(ChemicalProperties& CP, PhaseField& Phi, int x, int y, int z);
        std::vector<double> getWeightPercent(ChemicalProperties&CP, int x, int y, int z);
        void CalculateMoleFractions(ChemicalProperties& CP);
        void CalculateTotalMolarVolume(ChemicalProperties& CP);
        Tensor<double,1> AverageTotal();
        Tensor<double,1> AverageInterface(PhaseField& Phi);

        /*Tensor<double,3> AverageInterfacePair(PhaseField& Phi)
        {
            Tensor<double,3> result({Nphases,Nphases,Ncomp});
            Tensor<double,2> denom({Nphases,Nphases});

            result.set_to_zero();
            denom.set_to_zero();

            OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,Total,0,reduction(+:result,denom))
            if(Phi.Interface(x,y,z))
            {
                for(auto alpha = Phi.Fields(x,y,z).cbegin();
                        alpha != Phi.Fields(x,y,z).cend(); ++alpha)
                for(auto beta = Phi.Fields(x,y,z).cbegin();
                        beta != Phi.Fields(x,y,z).cend(); ++beta)
                if(alpha != beta)
                {
                    int locn = Phi.FieldsStatistics[alpha->index].Phase;
                    int locm = Phi.FieldsStatistics[beta->index].Phase;

                    for(int comp = 0; comp < Ncomp; comp++)
                    {
                        result({locn,locm,comp}) += Total(x,y,z)({comp});
                    }
                    denom({locn,locm}) += 1.0;
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END

            for(int n = 0; n < Nphases; n++)
            for(int m = 0; m < Nphases; m++)
            if(n != m)
            for(int comp = 0; comp < Ncomp; comp++)
            result({n,m,comp}) /= denom({n,m});

            return result;
        }*/

        Tensor<double,1> AveragePhase(PhaseField& Phi, int alpha);
        double get_CF(int phase, int constituent, int x, int y, int z);
        Tensor<double,1> get_CFs(int phase, int x, int y, int z);
        void set_CF(int phase, int constituent, double value, int x, int y, int z);
        void set_CFs(int phase, Tensor<double,1> value, int x, int y, int z);
        double get_MF(int phase, int component, int x, int y, int z);
        Tensor<double,1> get_MFs(int phase, int x, int y, int z);
        Tensor<double,2> get_MFs(int x, int y, int z);
        void set_MF(int phase, int component, double value, int x, int y, int z);
        void set_MFs(int phase, Tensor<double,1> value, int x, int y, int z);
    protected:
        int AtStart;
    private:
    };

    inline Tensor<double,1> Composition::AverageTotal()
    {
        Tensor<double,1> result({Ncomp});
        for(int comp = 0; comp < Ncomp; comp++)
        {
            double temp = 0.0;
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,Total,0,reduction(+:temp))
            {
                temp += Total(x,y,z)({comp});
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            result({comp}) = temp/(Nx*Ny*Nz);
        }

        return result;
    }

    inline Tensor<double,1> Composition::AverageInterface(PhaseField& Phi)
    {
        bool interfacepresent = false;
        Tensor<double,1> result({Ncomp});
        for(int comp = 0; comp < Ncomp; comp++)
        {
            int npoints = 0;
            double temp = 0.0;
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,Total,0,reduction(+:temp,npoints))
            if(Phi.Interface(x,y,z))
            {
                temp += Total(x,y,z)({comp});
                npoints++;
                interfacepresent = true;
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            result({comp}) = temp/(double)npoints;
        }

        if(!(interfacepresent))
        {
            result = AverageTotal();
        }

        return result;
    }

    inline Tensor<double,1> Composition::AveragePhase(PhaseField& Phi, int alpha)
    {
        Tensor<double,1> result({Ncomp});
        for(int comp = 0; comp < Ncomp; comp++)
        {
            double temp = 0.0;
            OMP_PARALLEL_STORAGE_LOOP_BEGIN(x,y,z,Total,0,reduction(+:temp))
            {
                double PhaseFractions = Phi.Fractions(x,y,z)({alpha});
                if(PhaseFractions > DBL_EPSILON)
                {
                    temp += Total(x,y,z)({comp})/PhaseFractions;
                }
            }
            OMP_PARALLEL_STORAGE_LOOP_END
            result({comp}) = temp;
        }

        return result;
    }

    inline double Composition::get_CF(int phase, int constituent, int x, int y, int z)
    {
        /** For quick access to the Constituent Fractions. */
        return ConstituentFractions[phase](x,y,z)({constituent});
    }

    inline Tensor<double,1> Composition::get_CFs(int phase, int x, int y, int z)
    {
        /** For quick access to all Constituent Fractions. */
        return ConstituentFractions[phase](x,y,z);
    }

    inline void Composition::set_CF(int phase, int constituent, double value, int x, int y, int z)
    {
        /** For quick write access to the Constituent Fractions.*/
        ConstituentFractions[phase](x,y,z)({constituent}) = value;
    }

    inline void Composition::set_CFs(int phase, Tensor<double,1> value, int x, int y, int z)
    {
        /** For quick write access to all Constituent Fractions.*/
        unsigned int Ncomp = value.size();
        if(Ncomp != ConstituentFractions[phase](x,y,z).size())
        {
            std::cout << "Mismatching Sitefraction storages sizes!" << std::endl;
            exit(1);
        }
        ConstituentFractions[phase](x,y,z) = value;
    }

    inline double Composition::get_MF(int phase, int component, int x, int y, int z)
    {
        /** For quick access to the Mole Fractions of a single phase. */
        return Phase(x,y,z)({phase,component});
    }

    inline Tensor<double,1> Composition::get_MFs(int phase, int x, int y, int z)
    {
        /** For quick access to all Mole Fractions of a single phase. */
        int Ncomp = Phase(x,y,z).size(1);
        Tensor<double,1> temp({Ncomp});
        for(int comp = 0; comp < Ncomp; comp++)
        temp({comp}) = Phase(x,y,z)({phase,comp});
        return temp;
    }

    inline Tensor<double,2> Composition::get_MFs(int x, int y, int z)
    {
        /** For quick access to all Mole Fractions of all phases. */
        return Phase(x,y,z);
    }

    inline void Composition::set_MF(int phase, int component, double value, int x, int y, int z)
    {
        /** For quick write access to the Mole Fractions.*/
        Phase(x,y,z)({phase,component}) = value;
    }

    inline void Composition::set_MFs(int phase, Tensor<double,1> value, int x, int y, int z)
    {
        /** For quick write access to all Mole Fractions of a single phase.*/
        unsigned int Ncomp = value.size();

        if(Ncomp != Phase(x,y,z).size(1))
        {
            exit(1);
        }
        for(int comp = 0; comp < (int)Ncomp; comp++)
        Phase(x,y,z)({phase,comp}) = value({comp});
    }
}
#endif