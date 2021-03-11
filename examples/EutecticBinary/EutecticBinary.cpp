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
    MySettings.Initialize();
    MySettings.ReadInput();
    Mb.Initialize(MySettings);
    Sigma.Initialize(MySettings);
    Phase.Initialize(MySettings);
    dG.Initialize(MySettings);
    Dob.Initialize(MySettings);
    Df.Initialize(MySettings);
    Tx.Initialize(MySettings);
    Bc.Initialize(MySettings);
    Cx.Initialize(MySettings);
    Tm.Initialize(MySettings, "去除时间统计");

    Df.ReadInput();
    Cx.ReadInput();
    dG.ReadInput();
    Tx.ReadInput();
    Bc.ReadInput();
    Mb.ReadInput();
    
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
        Initializations::Fractional(Phase, 0, 1, iR/2, Bc, MySettings);                //二元共晶体系的初始化
        Initializations::Sphere(Phase, 2, iR, (MySettings.Nx)/2, (MySettings.Ny)/2, 0, Bc, MySettings);           //圆形初晶相的初始化

        Df.SetInitialComposition(Phase, Cx);                    //设置初始状态的成分值
        Tx.SetInitial(Bc, Phase, 0);                            //设置初始状态的温度值
    }
    cout << "进入时间循环！" << endl;
    //-------------------------------进入循环----------------------------------//
    for(int tStart = MySettings.tStart + 1; tStart < MySettings.nSteps + 1; tStart++){

        Tm.SetStart();
            Sigma.CalculateCubic(Phase);
            Mb.CalculateCubic(Phase);
        Tm.SetTimeStamp("界面的能量和移动性计算");
            Df.SetPhaseFractions(Phase);
            double I_En = 0.0;
        Tm.SetTimeStamp("设置相分数");
            if(!(tStart%MySettings.tScreenWrite))
                I_En = Dob.Energy(Phase, Sigma);
            Dob.CalculatePhaseFieldIncrements(Phase, Sigma, Mb);
        Tm.SetTimeStamp("计算相场增量");
            Df.GetDrivingForce(Phase, Cx, Tx, dG);
        Tm.SetTimeStamp("得到驱动力");
            dG.Average(Phase, Bc);
        Tm.SetTimeStamp("驱动力平均");
            if(!(tStart%MySettings.tFileWrite))
            dG.WriteVTKforPhases(Phase, tStart);
            dG.MergePhaseFieldIncrements(Phase, Sigma, Mb);
        Tm.SetTimeStamp("合并相场增量");
            Phase.NormalizeIncrements(Bc, MySettings.dt);
        Tm.SetTimeStamp("增量标准化");
            Df.Solve(Phase, Cx, Tx, Bc, MySettings.dt);
        Tm.SetTimeStamp("求解计算");
            Tx.Set(Bc, MySettings.dt);
        Tm.SetTimeStamp("温度设置");
            Phase.MergeIncrements(Bc, MySettings.dt);
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
                     << "=============================\n" << endl;  
           
                dG.PrintDiagnostics();
                Phase.PrintPFVolumes();
                Tm.PrintWallClockSummary();                      //将时间步、界面能、自由能、体积等参数输出到屏幕上
            }
    } //结束循环
    return 0;
}