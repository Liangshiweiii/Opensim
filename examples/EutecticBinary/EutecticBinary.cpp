#include "Settings.h"
#include "Initializations.h"
#include "BoundaryConditions.h"
#include "Compositions.h"
#include "Diffusion.h"
#include "DoubleObstacle.h"
#include "InterfaceEnergy.h"
#include "InterfaceMobility.h"
#include "DrivingForce.h"
#include "PhaseField.h"
#include "Temperatures.h"
#include "TextOutput.h"
#include "Tools/TimeInfo.h"

using namespace std;
using namespace opensim;

int main(int argc, char *argv[]){
    //define modules
    Settings MySettings;
    Initializations Init;
    PhaseField Phase;
    BoundaryConditions Bc;
    TextOutput Out;//todo
    Temperature Tx;
    Composition Cx;
    InterfaceEnergy Sigma;
    InterfaceMobility Mb;
    DoubleObstacle Dob;
    DrivingForce dG;
    Diffusion Df;
    TimeInfo Tm;//todo
    //Initialize modules obj's initialization
    MySettings.Initialize();                  //基础设定模块的初始化
    MySettings.ReadInput();                   //基础设定模块从文件输入
    Mb.Initialize(MySettings);                //界面模块的初始化  
    Sigma.Initialize(MySettings);             //界面能模块的初始化 
    Phase.Initialize(MySettings);             //相场模块的初始化
    dG.Initialize(MySettings);                //驱动力模块的初始化
    Dob.Initialize(MySettings);               //双井势模块的初始化 
    Df.Initialize(MySettings);                //扩散模块的初始化
    Tx.Initialize(MySettings);                //温度模块的初始化
    Bc.Initialize(MySettings);                //边界条件模块的初始化
    Cx.Initialize(MySettings);                //成分模块的初始化
    Tm.Initialize(MySettings, "去除时间统计"); //时间模块的初始化

    Df.ReadInput();                           //扩散数据的输入
    Cx.ReadInput();                           //成分数据的输入
    dG.ReadInput();                           //驱动力数据的输入
    Tx.ReadInput();                           //温度数据的输入
    Bc.ReadInput();                           //边界数据的输入
    Mb.ReadInput();                           //界面扩散数据的输入
    
    cout << "初始化过程结束！" << endl;
    
    if (MySettings.Restart){
        cout << "重启数据读入中！";
        Phase.Read(Bc, MySettings.tStart);
        Cx.Read(MySettings.Cp, Bc, MySettings.tStart);
        Tx.Read(MySettings.tStart);
        cout << "完成！" << endl;


    }
    else{
        MySettings.tStart = -1;
        int iR = (MySettings.Nx)/4;
        //single lamella
        Initializations::Fractional(Phase, 0, 1, iR/2, Bc, MySettings);                                           //二元共晶体系z方向的空间初始化
        Initializations::Sphere(Phase, 2, iR, (MySettings.Nx)/2, (MySettings.Ny)/2, 0, Bc, MySettings);           //二元共晶体系x，y平面的空间初始化，选择的参考点是底边的中点

        Df.SetInitialComposition(Phase, Cx);                       //设置体系初始状态的成分值
        Tx.SetInitial(Bc, Phase, 0);                               //设置体系初始状态的温度值
    }
    cout << "进入循环！" << endl;
    //-------------------------------进入循环----------------------------------//
    for(int tStart = MySettings.tStart + 1; tStart < MySettings.nSteps + 1; tStart++){

        Tm.SetStart();                                              //程序开始时间计时
            Sigma.CalculateCubic(Phase);                            //求出最大接界面能的循环
            Mb.CalculateCubic(Phase);                               //求出最大界面移动性参数
        Tm.SetTimeStamp("界面的能量和移动性计算");
            Df.SetPhaseFractions(Phase);                            //设置每个点的热力学相分数
            double I_En = 0.0;
        Tm.SetTimeStamp("设置相分数");
            if(!(tStart%MySettings.tScreenWrite))                   //输出到屏幕的时间步判断
                I_En = Dob.Energy(Phase, Sigma);
            Dob.CalculatePhaseFieldIncrements(Phase, Sigma, Mb);    //计算单位时间相场的增量，即偏微分dphi/dt
        Tm.SetTimeStamp("计算相场增量");
            Df.GetDrivingForce(Phase, Cx, Tx, dG);                  //每个点相的组成以及驱动力dG的计算
        Tm.SetTimeStamp("得到驱动力");
            dG.Average(Phase, Bc);                                   //施加边界条件。将垂直于界面的驱动力平均化，使得界面稳定
        Tm.SetTimeStamp("设置边界条件，计算平均");
            if(!(tStart%MySettings.tFileWrite))                      //判断是否等于输出时间步输出计算结果到文件
            dG.WriteVTKforPhases(Phase, tStart);                     //将相场信息输入到VTK文件中
            dG.MergePhaseFieldIncrements(Phase, Sigma, Mb);          //计算相场增量，同时考虑驱动力以及相的限制
        Tm.SetTimeStamp("合并相场增量");
            Phase.NormalizeIncrements(Bc, MySettings.dt);             //将相场变量的数值标准化限制在（0,1）之间
        Tm.SetTimeStamp("增量标准化");
            Df.Solve(Phase, Cx, Tx, Bc, MySettings.dt);               //在相场，组成，边界条件同时作用下，求解相成分
        Tm.SetTimeStamp("求解计算");
            Tx.Set(Bc, MySettings.dt);                                //温度的迭代求解
        Tm.SetTimeStamp("温度设置");
            Phase.MergeIncrements(Bc, MySettings.dt);                 //计算合并增量后的各项数值，包括边界条件，flag位置，梯度数值，体积和分数。
        Tm.SetTimeStamp("合并相场增量");
            if (Phase.Fields((MySettings.Nx)/4, (MySettings.Ny)/2, (MySettings.Nz)/2)[0] <= 0.4){
                //Phase.MoveFrame(0,0,1, Bc);
                //Tx.MoveFrame(0,0,1, Bc);
                //Cx.MoveFrame(0,0,1, Bc);
            }
            if (!(tStart%MySettings.tFileWrite)){
                //  输出VTK形式
                Phase.WriteVTK(tStart);
                Cx.WriteVTK(tStart);
                Cx.WriteStatistics(tStart, MySettings.dt);
            }
            if (!(tStart%MySettings.tRestartWrite)){
                //  输出到文本文件
                Phase.Write(tStart);
                Cx.Write(tStart);
                Tx.Write(tStart);
            }
        Tm.SetTimeStamp("输出到文件和屏幕");
            time_t rawtime;                                         
            time(&rawtime);

            if (!(tStart%MySettings.tScreenWrite)){
                cout << "++++++++++++++++++++++++++++\n"
                     << "时间步    : " << tStart << "\n"
                     << "实际时间  : " << ctime(&rawtime) << "\n";
                cout << "================================\n"               
                     << "界面能    : " << I_En << "\n"
                     << "=============================\n" << endl;  //将时间等计算信息输出到屏幕
           
                dG.PrintDiagnostics();                              //判断驱动力极限的函数
                Phase.PrintPFVolumes();                             //输出相的空间体积
                Tm.PrintWallClockSummary();                         //将时间步、界面能、自由能、体积等参数输出到屏幕上
            }
    } //结束循环
    return 0;
}