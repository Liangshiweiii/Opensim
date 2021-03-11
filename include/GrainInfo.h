#ifndef GRAININFO_H
#define GRAININFO_H

#include "Tools/Includes.h"
#include "Tools/UserInterface.h"
#include "Tools/Quaternion.h"
#include "Settings.h"
namespace opensim
{

class Grain                                                                     ///<存储关于grain (phasefield)的信息，用作GrainInfo类中的存储空间
{
 public:
    long double MAXVolume;                                                      ///< 模拟过程中相场的最大体积
    long double Volume;                                                         ///< 空间的总容积
    double      Density;                                                        ///< 晶粒平均密度
    double      WettingParameter;                                               ///< LBM中本兹力的润湿参数
    dVector3    Rcm;                                                            ///< 质心的坐标
    dVector3    Vcm;                                                            ///< 质心的速度
    dVector3    Acm;                                                            ///< 质心的加速度
    dVector3    aVel;                                                           ///< 角速度
    dVector3    aAcc;                                                           ///< 角加速度
    dVector3    Force;                                                          ///< 作用在纹理上的力
    dVector3    Torque;                                                         ///< 施加在晶粒上的扭矩
    dMatrix3x3  InertiaM;                                                       ///< 转动惯量张量
    Quaternion  Orientation;                                                    ///< 纹理相对于参照系的方向

    int Exist;                                                                  ///< 如果仿真域中存在相场，则为“真”
    int IsMobile;                                                               ///< 如果晶粒能移动，那就对了。与固液耦合运动相关
    int Phase;                                                                  ///< 给定相场的热力学相
    int Stage;                                                                  ///< 晶粒阶段:1 ->核，0 ->全尺寸晶粒。
    int State;                                                                  ///< 相(固体、液体或气体)的聚集状态。
    int Variant;                                                                ///< 给定热力学相的对称变体

    void Read(std::fstream& In)                                                 ///< 从给定的文件流读取晶体信息
    {
        In >> State;
        In >> Stage;
        In >> Exist;
        In >> IsMobile;
        In >> Phase;
        In >> Variant;
        In >> Volume;
        In >> Density;
        In >> Rcm[0] >> Rcm[1] >> Rcm[2];
        In >> Vcm[0] >> Vcm[1] >> Vcm[2];
        In >> Acm[0] >> Acm[1] >> Acm[2];
        In >> aVel[0] >> aVel[1] >> aVel[2];
        In >> aAcc[0] >> aAcc[1] >> aAcc[2];
        In >> Force[0] >> Force[1] >> Force[2];
        In >> Torque[0] >> Torque[1] >> Torque[2];
        double Q[4];
        In >> Q[0] >> Q[1] >> Q[2] >> Q[3];
        Orientation.set(Q[0], Q[1], Q[2], Q[3]);
    }

    void Write(std::ostream& Out) const                                         ///< 将晶体信息写入给定的文件流
    {
        Out << State    << std::endl;
        Out << Stage    << std::endl;
        Out << Exist    << std::endl;
        Out << IsMobile << std::endl;
        Out << Phase    << std::endl;
        Out << Variant  << std::endl;
        Out << Volume   << std::endl;
        Out << Density  << std::endl;
        Out << Rcm.getX()       << " " << Rcm.getY()       << " " << Rcm.getZ()       << std::endl;
        Out << Vcm.getX()       << " " << Vcm.getY()       << " " << Vcm.getZ()       << std::endl;
        Out << Acm.getX()       << " " << Acm.getY()       << " " << Acm.getZ()       << std::endl;
        Out << aVel.getX()      << " " << aVel.getY()      << " " << aVel.getZ()      << std::endl;
        Out << aAcc.getX()      << " " << aAcc.getY()      << " " << aAcc.getZ()      << std::endl;
        Out << Force.getX()     << " " << Force.getY()     << " " << Force.getZ()     << std::endl;
        Out << Torque.getX()    << " " << Torque.getY()    << " " << Torque.getZ()    << std::endl;
        Out << Orientation[0]   << " " << Orientation[1]   << " " << Orientation[2]   << " " << Orientation[3] << std::endl;
    }

    void WriteTable(std::ostream& Out, const int tStep, const char sep = ',') const ///< 将信息写入csv输出
    {
        // 写标题
        if (tStep == 0)
        {
            Out  << "Time"       << sep
                 << "Volume"     << sep
                 << "Density"    << sep
                 << "X"          << sep << "Y"         << sep << "Z"         << sep
                 << "VX"         << sep << "VX"        << sep << "VZ"        << sep
                 << "AngularVX"  << sep << "AngularVY" << sep << "AngularVZ" << sep
                 << "AX"         << sep << "AX"        << sep << "AZ"        << sep
                 << "AngularAZ"  << sep << "AngularAY" << sep << "AngularAZ" << sep
                 << "ForceX"     << sep << "ForceY"    << sep << "ForceZ"    << sep
                 << "TorqueX"    << sep << "TorqueY"   << sep << "TorqueZ"   << sep
                 << "Q1"         << sep << "Q2"        << sep << "Q3"        << sep << "Q4"<< std::endl;
        }

        Out << tStep            << sep;
        Out << Volume           << sep;
        Out << Density          << sep;
        Out << Rcm.getX()       << sep << Rcm.getY()       << sep << Rcm.getZ()       << sep;
        Out << Vcm.getX()       << sep << Vcm.getY()       << sep << Vcm.getZ()       << sep;
        Out << aVel.getX()      << sep << aVel.getY()      << sep << aVel.getZ()      << sep;
        Out << Acm.getX()       << sep << Acm.getY()       << sep << Acm.getZ()       << sep;
        Out << aAcc.getX()      << sep << aAcc.getY()      << sep << aAcc.getZ()      << sep;
        Out << Force.getX()     << sep << Force.getY()     << sep << Force.getZ()     << sep;
        Out << Torque.getX()    << sep << Torque.getY()    << sep << Torque.getZ()    << sep;
        Out << Orientation[0]   << sep << Orientation[1]   << sep << Orientation[2]   << sep << Orientation[3] << std::endl;
    }

    Grain()                                                                     ///< 构造函数-设置所有值为默认值
    {
        Clear();
    }

    void Clear()                                                                ///< 清除晶粒信息
    {
        State = Solid;
        Stage = 0;
        Exist = 0;
        IsMobile = 0;
        Phase = 0;
        Variant = 0;
        Volume  = 0.0;
        MAXVolume  = 0.0;
        Density = 0.0;

        Orientation.set_to_zero();

        InertiaM.set_to_zero();
        Rcm.set_to_zero();
        Vcm.set_to_zero();
        Acm.set_to_zero();
        aVel.set_to_zero();
        aAcc.set_to_zero();
        Torque.set_to_zero();
        Force.set_to_zero();
    };
    bool IsNucleus(void)
    {
        return Stage;
    };
    bool IsPresent(void)
    {
        return Exist;
    };
};

class GrainInfo
{
 public:

    void Allocate(const int size)
    {
        GrainStorage.resize(size);
        for(unsigned int idx = 0; idx < GrainStorage.size(); idx++)
        {
            GrainStorage[idx].Clear();
        }
    }
    Grain& operator[](const int index)
    {
#ifdef DEBUG
        if(index >= static_cast<int>(GrainStorage.size()))
        {
            std::cerr << "Error in GrainInfo::operator[]\n"
                      << "Access beyond storage range. index = "
                      << index << " >= size = "
                      << GrainStorage.size()
                      << "\nTerminating!!!" << std::endl;
            std::cout.flush();
            exit(13);
        }
#endif
        return GrainStorage[index];
    }
    Grain const& operator[](const int index) const
    {
#ifdef DEBUG
        if(index >= static_cast<int>(GrainStorage.size()))
        {
            std::cerr << "Error in GrainInfo::operator[]\n"
                      << "Access beyond storage range. index = "
                      << index << " >= size = "
                      << GrainStorage.size()
                      << "\nTerminating!!!" << std::endl;
            std::cout.flush();
            exit(13);
        }
#endif
        return GrainStorage[index];
    }

    size_t size() const
    {
        return GrainStorage.size();
    }

    int add_grain(const int PhaseIndex)
    {
        Grain locGrain;
        locGrain.Clear();
        locGrain.Exist = true;
        locGrain.Phase = PhaseIndex;
        GrainStorage.push_back(locGrain);
        return GrainStorage.size() - 1;
    }

    int add_nucleus(const int PhaseIndex)
    {
        for(unsigned int i = 0; i != GrainStorage.size(); i++)
        if(not GrainStorage[i].Exist)
        {
            GrainStorage[i].Exist = true;
            GrainStorage[i].Stage = 1;
            GrainStorage[i].Phase = PhaseIndex;
            return i;
        }
        Grain locGrain;
        locGrain.Clear();
        locGrain.Exist = true;
        locGrain.Stage = 1;
        locGrain.Phase = PhaseIndex;
        GrainStorage.push_back(locGrain);

        return GrainStorage.size() - 1;
    }
    void ChangePhaseIndex(int PhaseOld, int PhaseNew)
    {
        for(unsigned int idx = 0; idx < GrainStorage.size(); idx++)
        if(GrainStorage[idx].Phase == PhaseOld)
        {
            GrainStorage[idx].Phase = PhaseNew;
        }
    }
    void ChangePhaseIndices(std::vector<int> PhaseOld, std::vector<int> PhaseNew)
    {
        if(PhaseOld.size() != PhaseNew.size())
        {
            std::cout << "GrainInfo::ChangePhase(initilizer_list)!" << std::endl
                      << "Incompatible phase indeces list sizes: "
                      << PhaseOld.size() << "!=" << PhaseNew.size() << "!" << std::endl;
            exit(1);
        }
        for(unsigned int idx = 0; idx < GrainStorage.size(); idx++)
        for(unsigned int n = 0; n < PhaseOld.size();n++)
        if(GrainStorage[idx].Phase == PhaseOld[n])
        {
            GrainStorage[idx].Phase = PhaseNew[n];
        }
    }
    void Write(const int tStep) const
    {
        std::string FileName = UserInterface::MakeFileName(RawDataDir, "GrainStatistics_", tStep,".dat");

        std::fstream GrainInfoOut(FileName, std::ios::out);
        GrainInfoOut << GrainStorage.size() << std::endl;
        for(unsigned int idx = 0; idx < GrainStorage.size(); idx++)
        {
            GrainStorage[idx].Write(GrainInfoOut);
        }
        GrainInfoOut.close();
    }

    void Read(const std::string FileName)
    {
        std::fstream GrainInfoIn(FileName, std::ios::in);

        unsigned int size = 0;
        GrainInfoIn >> size;
        if(GrainStorage.size() < size)
        {
            GrainStorage.resize(size);
        }

        for(unsigned int idx = 0; idx < GrainStorage.size(); idx++)
        {
            GrainStorage[idx].Read(GrainInfoIn);
        }
        GrainInfoIn.close();
    }

    void Read(const int tStep)
    {
        std::string FileName = UserInterface::MakeFileName(RawDataDir, "GrainStatistics_", tStep,".dat");

        Read(FileName);
    }

    void Write(int ID, int tStep)
    {
        std::ostringstream FileName;
        FileName.fill('0');
        FileName << "Checkpoint/PhaseField/GrainStats" << std::setw(8) << ID << "_" << tStep << ".dat";

        std::fstream GrainInfoOut(FileName.str(), std::ios::out);
        GrainInfoOut << GrainStorage.size() << std::endl;
        for(unsigned int idx = 0; idx < GrainStorage.size(); idx++)
        {
            GrainStorage[idx].Write(GrainInfoOut);
        }
        GrainInfoOut.close();
    }

    void Read(int ID, int tStep)
    {
        std::ostringstream FileName;
        FileName.fill('0');
        FileName << "Checkpoint/PhaseField/GrainStats" << std::setw(8) << ID << "_" << tStep << ".dat";

        std::fstream GrainInfoIn(FileName.str(), std::ios::in);

        int size = 0;
        GrainInfoIn >> size;
        GrainStorage.resize(size);

        for(unsigned int idx = 0; idx < GrainStorage.size(); idx++)
        {
            GrainStorage[idx].Read(GrainInfoIn);
        }
        GrainInfoIn.close();
    }
    std::vector < Grain > GrainStorage;
 protected:
 private:
};

}// namespace openphase
#endif
