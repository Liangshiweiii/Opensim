
#include "Info.h"
#include "DrivingForce.h"
#include "Tools/UserInterface.h"
#include "VTK.h"
#include "GrainInfo.h"
#include "Settings.h"
#include "PhaseField.h"
#include "InterfaceMobility.h"
#include "InterfaceEnergy.h"
#include "BoundaryConditions.h"
#include <random>
namespace opensim
{

using namespace std;
DrivingForce::DrivingForce(const Settings& locSettings, const int boundary)
{
    this->Initialize(locSettings, boundary);
    this->ReadInput();
}

void DrivingForce::Initialize(const Settings& locSettings,
        const int boundary)
{
    thisclassname = "DrivingForce";
    //DefaultInputFileName = ProjectInputDir + "DrivingForceInput.opi";

    Nx = locSettings.Nx;
    Ny = locSettings.Ny;
    Nz = locSettings.Nz;

    // Setting default values to be used if ReadInput() is not called
    CutOff = 0.95;
    Averaging = false;
    Range = ((locSettings.iWidth + 1)/2)+1;
    PhiThreshold = 1.0/3.0;
    switch(int(locSettings.iWidth))
    {
        case 3:
        {
            PhiThreshold = 1.0/5.0;
            break;
        }
        case 4:
        {
            PhiThreshold = 1.0/4.0;
            break;
        }
        case 5:
        {
            PhiThreshold = 1.0/3.0;
            break;
        }
        default: PhiThreshold = 1.0/3.0;
    }
    // End setting default values

    Raw.Allocate(Nx, Ny, Nz, Range);

    dGabLimitReachedCounter = 0;
    EnDensUnits = locSettings.UnitsOfEnergy + "/" + locSettings.UnitsOfLength + "^3";

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}
void DrivingForce::ReadInput(string InputFileName)
{
    fstream inp(InputFileName.c_str(), ios::in);

    if (!inp)
    {
        Info::WriteExit("File \"" + InputFileName + "\" could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };

    Info::WriteBlankLine();

    Info::WriteLineInsert("Driving Force settings");
    Info::WriteStandard("Source", InputFileName.c_str());

    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);

    CutOff = UserInterface::ReadParameterD(inp, moduleLocation, string("CutOff"), false, CutOff);           //系数，局部多少的驱动力被去除以维持界面的稳定
    Averaging = UserInterface::ReadParameterB(inp, moduleLocation, string("Averg"), false, "Yes");
    Range = UserInterface::ReadParameterI(inp, moduleLocation, string("Range"), false, Range);              
    PhiThreshold = UserInterface::ReadParameterD(inp, moduleLocation, string("Thresh"), false, PhiThreshold);

    Info::WriteLine();

    inp.close();
}

void DrivingForce::Clear()
{
    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,Raw.Bcells(),)
    {
        Raw(i,j,k).clear();
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DrivingForce::SetBoundaryConditions(BoundaryConditions& BC)
{
    /*if(Nx > 1)*/ BC.SetX(Raw);
    /*if(Ny > 1)*/ BC.SetY(Raw);
    /*if(Nz > 1)*/ BC.SetZ(Raw);
}

void DrivingForce::Average(PhaseField& Phase, BoundaryConditions& BC)
{
    if(Averaging)
    {
        SetBoundaryConditions(BC);
        CollectAverage(Phase);
        SetBoundaryConditions(BC);
        DistributeAverage(Phase);
    }
}

void DrivingForce::CollectAverage(PhaseField& Phase)
{
    const int Xrange = min(Range, Nx-1);
    const int Yrange = min(Range, Ny-1);
    const int Zrange = min(Range, Nz-1);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        if (Phase.Interface(i,j,k))
        for(auto it = Raw(i,j,k).begin();
                 it < Raw(i,j,k).end(); ++it)
        {
            double PhiAlpha = Phase.Fields(i,j,k).get(it->indexA);
            double PhiBeta  = Phase.Fields(i,j,k).get(it->indexB);

			if (PhiBeta > 0 && PhiAlpha > 0)
            if (PhiAlpha/PhiBeta > PhiThreshold/(1.0 - PhiThreshold) and
                PhiAlpha/PhiBeta < (1.0 - PhiThreshold)/PhiThreshold)
            {
                double value = 0.0;

                double SumWeights = 0.0;

                for(int ii = -Xrange; ii <= Xrange; ii++)
                for(int jj = -Yrange; jj <= Yrange; jj++)
                for(int kk = -Zrange; kk <= Zrange; kk++)
                if(Phase.Interface(i+ii,j+jj,k+kk))                     //确认点是否在界面上
                {
                    double dist = sqrt(ii*ii + jj*jj + kk*kk);          //计算到中心点的距离distance

                    double locPhiAlphaValue = Phase.Fields(i+ii, j+jj, k+kk)[it->indexA];
                    double locPhiBetaValue = Phase.Fields(i+ii, j+jj, k+kk)[it->indexB];

                    double weight = sqrt(locPhiAlphaValue * locPhiBetaValue)*(Range - dist);
                    if (weight > DBL_EPSILON)
                    {
                        SumWeights += weight;
                        value += weight*Raw(i+ii, j+jj, k+kk).get_asym(it->indexA, it->indexB);
                    }
                }
                if(SumWeights > 0.0) it->value2 = value/SumWeights;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DrivingForce::DistributeAverage(PhaseField& Phase)
{
    const int Xrange = min(Range, Nx-1);
    const int Yrange = min(Range, Ny-1);
    const int Zrange = min(Range, Nz-1);

    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,0,)
    {
        if (Phase.Interface(i,j,k))
        for(auto it = Raw(i,j,k).begin();
                 it < Raw(i,j,k).end(); ++it)
        {
            double PhiAlpha = Phase.Fields(i,j,k).get(it->indexA);
            double PhiBeta  = Phase.Fields(i,j,k).get(it->indexB);

            if (PhiAlpha*PhiBeta != 0.0)
            {
                double counter = 0.0;
                double value = 0.0;
                for(int ii = -Xrange; ii <= Xrange; ii++)
                for(int jj = -Yrange; jj <= Yrange; jj++)
                for(int kk = -Zrange; kk <= Zrange; kk++)
                if(Phase.Interface(i+ii,j+jj,k+kk))
                {
                    double dist = sqrt(ii*ii + jj*jj + kk*kk);
                    double locPhiAlpha = Phase.Fields(i+ii,j+jj,k+kk).get(it->indexA);
                    double locPhiBeta  = Phase.Fields(i+ii,j+jj,k+kk).get(it->indexB);
					if (locPhiBeta > 0 && locPhiAlpha > 0)
                    if (locPhiAlpha/locPhiBeta > PhiThreshold/(1.0 - PhiThreshold) and
                        locPhiAlpha/locPhiBeta < (1.0 - PhiThreshold)/PhiThreshold and (Range - dist) > 0.0)
                    {
                        counter ++;
                        value += Raw(i+ii, j+jj, k+kk).get_asym2(it->indexA, it->indexB);
                    }
                }
                if(counter > 0.0)
                {
                    it->value = value/counter;
                }
            }
            else if((!Phase.FieldsStatistics[it->indexA].Stage and
                     !(Phase.FieldsStatistics[it->indexA].MAXVolume == 0.0)) and
                    (!Phase.FieldsStatistics[it->indexB].Stage and
                     !(Phase.FieldsStatistics[it->indexB].MAXVolume == 0.0)))
            {
                it->value = 0.0;
            }
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
}

void DrivingForce::MergePhaseFieldIncrements(PhaseField& Phase,
                                             InterfaceEnergy& Sigma,
                                             InterfaceMobility& Mu)
{
    if(CutOff == 0.0)                 //如果驱动力为0时重新读取输入参数
    {
        stringstream message;
        message << "CutOff = 0.0, please correct the driving force input file "
                << "or add missing call to DrivingForce::ReadInput() in your code!";
        Info::WriteExit(message.str(), thisclassname, "MergePhaseFieldIncrements()");          //输出终止时的函数名
        exit(13);
    }

    const double Eta = Phase.Eta;
    double locMaxPsi = 0.0;
	double locMAXdGabOvershoot = MAXdGabOvershoot;
	double locMAXDrivingForce = MAXDrivingForce;
	int locDGabLimitReachedCounter = dGabLimitReachedCounter;

	OMP_PARALLEL_STORAGE_LOOP_BEGIN(i, j, k, Phase.Fields, 0, reduction(+:locDGabLimitReachedCounter) reduction(max:locMAXdGabOvershoot, locMAXDrivingForce, locMaxPsi)) // reduction(+:dGabLimitReachedCounter)
    {
        if (Phase.Interface(i,j,k))
        {
            double Prefactor = 2.0*Pi/(Eta*double(Phase.Fields(i,j,k).size()));
            for(auto it = Raw(i,j,k).cbegin();
                     it < Raw(i,j,k).cend(); ++it)
            {
                double alphaVal = Phase.Fields(i,j,k).get(it->indexA);
                double  betaVal = Phase.Fields(i,j,k).get(it->indexB);

                double norm = sqrt(alphaVal * betaVal);

                if((norm <= DBL_EPSILON) and
                   ((Phase.FieldsStatistics[it->indexA].Stage and Phase.FieldsStatistics[it->indexA].MAXVolume == 0.0) or
                    (Phase.FieldsStatistics[it->indexB].Stage and Phase.FieldsStatistics[it->indexB].MAXVolume == 0.0)))
                {
                    norm = 1.0e-6;
                }

                /*if(Phase.FieldsStatistics[it->indexA].Stage)
                {
                    norm += 0.5*(1.0 - Phase.FieldsStatistics[it->indexA].Volume/Phase.RefVolume);
                }
                if(Phase.FieldsStatistics[it->indexB].Stage)
                {
                    norm += 0.5*(1.0 - Phase.FieldsStatistics[it->indexB].Volume/Phase.RefVolume);
                }*/

                double locCutOff = CutOff;
                if(Phase.FieldsStatistics[it->indexA].Stage or Phase.FieldsStatistics[it->indexB].Stage)
                {
                    locCutOff += 1.0 - min(Phase.FieldsStatistics[it->indexA].Volume,
                                           Phase.FieldsStatistics[it->indexB].Volume)/Phase.RefVolume;
                }

                double AllowedDrivingForce = locCutOff*Prefactor*
                                Sigma.MinIntEnergy(Phase.FieldsStatistics[it->indexA].Phase,
                                                   Phase.FieldsStatistics[it->indexB].Phase);

                double dPsi_dt  = AllowedDrivingForce *
                                  tanh( it->value  / AllowedDrivingForce); //dGmax*Tanh(dG/dGmax)

                dPsi_dt *= Mu(i,j,k, it->indexA, it->indexB) * norm * Prefactor;

				locMaxPsi = max(fabs(dPsi_dt), locMaxPsi);

                Phase.FieldsDot(i,j,k).add_asym(it->indexA, it->indexB,  dPsi_dt);
                // Begin collecting statistics:
                double dGloc = fabs(it->value);
                if(dGloc > 0.4*AllowedDrivingForce)
                {
					locDGabLimitReachedCounter++;
                }

                double tmpMAXdGabOvershoot  = dGloc/AllowedDrivingForce;
                if(tmpMAXdGabOvershoot > locMAXdGabOvershoot)
                {
					locMAXdGabOvershoot = tmpMAXdGabOvershoot;
					locMAXDrivingForce = dGloc;
                }
                // End collecting statistics
            }
            Raw(i,j,k).clear();
        }
        else if(Phase.Fields(i,j,k).flag)
        {
            Raw(i,j,k).clear();
        }
    }
    OMP_PARALLEL_STORAGE_LOOP_END
	maxPsi = locMaxPsi;
	MAXdGabOvershoot = locMAXdGabOvershoot;
	MAXDrivingForce = locMAXDrivingForce;
	dGabLimitReachedCounter = locDGabLimitReachedCounter;
}

double DrivingForce::MaxTimeStep(PhaseField& Phase, InterfaceEnergy& Sigma,
                                InterfaceMobility& Mu, Settings& OPSettings,
                                double TheorLimit, double NumerLimit)
{
    /** This function will calculate the maximum time step the phase field
    should be solved with.
    This function uses two different stability criteria. A theoretical stability
    function as 0.25*dx*dx/mu/Sigma (when TheorLimit is set to 0.25), and a
    numerical criteria, where the maximum phase field increment in the whole
    simulation domain has to be lower than 1E-3 (if NumerLimit is set to 1E-3).
    If phase field gets unstable, decrease TheorLimit, if diffusion field
    becomes unstable due to too high phase field increments, decrease
    NumerLimit.*/

    double maxmu = Mu.maxMu;
    double maxsg = Sigma.maxSigma;
    double dx2 = OPSettings.dx*OPSettings.dx;
    double maxTheorTimeStep = 0.0;
    double maxNumerTimeStep = 0.0;

    // Calculate theoretical maximum allowed time step
    if((maxmu > DBL_EPSILON) and (maxsg > DBL_EPSILON))
    maxTheorTimeStep = TheorLimit*dx2/maxmu/maxsg;

    // Calculate maximum numerical allowed time step
    if(maxPsi > DBL_EPSILON)
    maxNumerTimeStep = NumerLimit/maxPsi;

    /* the equation would look like: maxNumerTimeStep=dt*NumerLimit/maxPsi
    if maxPsi would be the phase field increment. But as Psi it has to be
    multiplied with dt, it cancels out in this equation. So no dt needed!*/

    // Calculate total maximum allowed time step
    double maxTimeStep = min(maxTheorTimeStep,maxNumerTimeStep);

    return maxTimeStep;
}

void DrivingForce::AddNoiseRelative(double amplitude, int alpha, int beta)
{
    std::mt19937_64 NoiseRandGenerator(98310);
    std::uniform_real_distribution <double> NoiseDistribution(-amplitude, amplitude);
    for(int i = 0; i < Nx; ++i)
    for(int j = 0; j < Ny; ++j)
    for(int k = 0; k < Nz; ++k)
    for(auto it = Raw(i,j,k).begin(); it < Raw(i,j,k).end(); it++)
    if((it->indexA == alpha and it->indexB == beta) or 
       (it->indexB == alpha and it->indexA == beta))
    {
        it->value = it->value * (1.0 + NoiseDistribution(NoiseRandGenerator));
    }
}

void DrivingForce::AddNoiseAbsolute(double amplitude, int alpha, int beta)
{
    std::mt19937_64 NoiseRandGenerator(543454);
    std::uniform_real_distribution <double> NoiseDistribution(-amplitude, amplitude);
    for(int i = 0; i < Nx; ++i)
    for(int j = 0; j < Ny; ++j)
    for(int k = 0; k < Nz; ++k)
    for(auto it = Raw(i,j,k).begin(); it < Raw(i,j,k).end(); it++)
    if((it->indexA == alpha and it->indexB == beta) or 
       (it->indexB == alpha and it->indexA == beta))
    {
        it->value += NoiseDistribution(NoiseRandGenerator);
    }
}

void DrivingForce::PrintDiagnostics()
{
    if (dGabLimitReachedCounter)
    {
        std::string message = "\nThe driving force value has been limited " +
                               std::to_string(static_cast <long int> (dGabLimitReachedCounter)) + " times\n"
                               + "Max value: " + std::to_string(static_cast <long double> (MAXDrivingForce))
                               + " " + EnDensUnits + "\n"
                               + "Max overshoot: " + std::to_string(static_cast <long double> (MAXdGabOvershoot)) + " times\n";
        Info::WriteStandard(thisclassname + "::PrintDiagnostics()", message, false);
        dGabLimitReachedCounter = 0;
        MAXdGabOvershoot        = 0;
        MAXDrivingForce         = 0;
    }
    else
    {
        std::string message = "No Driving Force limiting was needed!";
        Info::WriteSimple(message, false);
    }
}

void DrivingForce::PrintPointStatistics(const int x, const int y, const int z) const
{
    std::stringstream pointstat;
    pointstat << "DrivingForce Indices:\t";

    for (auto alpha = Raw(x,y,z).cbegin();
              alpha < Raw(x,y,z).cend(); ++alpha)
    {
        pointstat << alpha->indexA << ", " << alpha->indexB << "\t\t";
    }
    pointstat << endl;
    pointstat << "DrivingForce Values:\t";

    for (auto alpha = Raw(x,y,z).cbegin();
              alpha < Raw(x,y,z).cend(); ++alpha)
    {
        pointstat << Raw(x,y,z).get_asym(alpha->indexA, alpha->indexB) << "\t\t";
    }
    pointstat << endl;

    Info::WriteSimple(pointstat.str(), false);
}

void DrivingForce::WriteVTK(const int tStep, const int indexA, const int indexB) const
{
    stringstream converter;
    converter << indexA << "_" << indexB << "_" ;
    string phases = converter.str();

    stringstream buffer;
    std::vector<int> DataTypes {PDScalars};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);
    {
        buffer << "<DataArray type = \"Float64\" Name = \"dG("
               << indexA << "," << indexB << ")\" "
               << "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            buffer << Raw(i,j,k).get_asym(indexA, indexB) << "\n";
        }
        buffer << "</DataArray>" << endl;
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "DrivingForce_" + phases, tStep);
}

void DrivingForce::WriteVTKforPhases(PhaseField& Phi, const int tStep) const
{
    /**This function will write a readable output in VTS format, summarizing
     * the driving forces of all individual grains as the driving force between
     * individual phase fields.*/

    stringstream converter;                 //stringstream是一个用于字符串输入输出的流
    string phases = converter.str();

    stringstream buffer;
    std::vector<int> DataTypes {PDScalars};

    VTK::WriteHeader(buffer, Nx, Ny, Nz);
    VTK::WriteBeginPointData(buffer, DataTypes);

    for(int alpha = 0; alpha < Phi.Nphases; alpha++)
    for(int beta = alpha+1; beta < Phi.Nphases; beta++)
    {
        buffer << "<DataArray type = \"Float64\" Name = \"dG("
               << alpha << "," << beta << ")\" "
               << "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < Nz; ++k)
        for(int j = 0; j < Ny; ++j)
        for(int i = 0; i < Nx; ++i)
        {
            double tempdG = 0.0;
            for(auto it1 = Phi.Fields(i,j,k).cbegin();
                     it1 < Phi.Fields(i,j,k).cend(); ++it1)
            for(auto it2 = Phi.Fields(i,j,k).cbegin();
                     it2 < Phi.Fields(i,j,k).cend(); ++it2)
            if((Phi.FieldsStatistics[it1->index].Phase == alpha)
            and(Phi.FieldsStatistics[it2->index].Phase == beta))
            {
                tempdG += Raw(i,j,k).get_asym(it1->index, it2->index);
            }
            buffer << tempdG << "\n";
        }
        buffer << "</DataArray>" << endl;
    }
    VTK::WriteEndPointData(buffer);
    VTK::WriteCoordinates(buffer, Nx, Ny, Nz);
    VTK::WriteToFile(buffer, "DrivingForce", tStep);               //将buffer中缓冲的字符串数据输入到文件中
}

void DrivingForce::Remesh(int newNx, int newNy, int newNz)
{
    Raw.Reallocate(newNx, newNy, newNz);

    Nx = newNx;
    Ny = newNy;
    Nz = newNz;

    Info::WriteStandard(thisclassname, "Remeshed");
}

DrivingForce& DrivingForce::operator= (const DrivingForce& rhs)
{
    // protect against invalid self-assignment and copy of unitialized object
    if (this != &rhs and rhs.thisclassname == "DrivingForce")
    {
        thisclassname = rhs.thisclassname;
        //DefaultInputFileName = rhs.DefaultInputFileName;
        Nx = rhs.Nx;
        Ny = rhs.Ny;
        Nz = rhs.Nz;
        Range = rhs.Range;
        PhiThreshold = rhs.PhiThreshold;

        dGabLimitReachedCounter = rhs.dGabLimitReachedCounter;
        MAXdGabOvershoot = rhs.MAXdGabOvershoot;
        MAXDrivingForce = rhs.MAXDrivingForce;
        EnDensUnits = rhs.EnDensUnits;
        CutOff = rhs.CutOff;
        Averaging = rhs.Averaging;
        EnDensUnits = rhs.EnDensUnits;

        if (Raw.IsNotAllocated())
        {
            Raw.Allocate(Nx, Ny, Nz, rhs.Raw.Bcells());
        }
        else if (not Raw.IsSize(rhs.Nx, rhs.Ny, rhs.Nz))
        {
            Raw.Reallocate(Nx, Ny, Nz);
        }

        OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Raw,Raw.Bcells(),)
        {
            Raw(i,j,k) = rhs.Raw(i,j,k);
        }
        OMP_PARALLEL_STORAGE_LOOP_END
    }
    return *this;
}

}// namespace opensim


