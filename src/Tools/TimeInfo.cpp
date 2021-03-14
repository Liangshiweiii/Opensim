#include "Tools/TimeInfo.h"
#include "Info.h"

namespace opensim
{
    using namespace std;

    TimeInfo::TimeInfo(const Settings& locSettings, const std::string Name)
    {
        this->Initialize(locSettings, Name);
    }

    TimeInfo::~TimeInfo(void)
    {
        Info::WriteStandard(thisclassname, "Exited normally");
    }

    void TimeInfo::Initialize(const Settings& locSettings, const string Name)
    {
        thisclassname = "TimeInfo";
        TimerName = Name;
    }
    void TimeInfo::Reset(void)
    {
        TimeCollector.clear();
        Messages.clear();
    }

    void TimeInfo::SetStart(void)
    {
        Reset();
        TimeCollector.push_back(clock());
        Messages.push_back("Start");
    }

    void TimeInfo::SetTimeStamp(const string Message)
    {
        TimeCollector.push_back(clock());
        Messages.push_back(Message);
    }

    void TimeInfo::PrintWallClockSummary(void) const
    {
        Info::WriteLine("=");
        Info::WriteSimple(TimerName);
        Info::WriteLine("-");

        clock_t TotalTime = 0;
        for(unsigned int it = 1; it < TimeCollector.size(); it++)
        {
            if(Messages[it]!="#")
            {
                TotalTime += (TimeCollector[it] - TimeCollector[it-1]);
            }
        }
        double TotalConsumedTime = double(TotalTime) / double(CLOCKS_PER_SEC*omp_get_max_threads());

        for(unsigned int it = 1; it < TimeCollector.size(); it++)
        {
            if(Messages[it]!="#")
            {
                stringstream showtime;
                double SectionConsumedTime = double(TimeCollector[it] - TimeCollector[it-1]) / double(CLOCKS_PER_SEC*omp_get_max_threads());
                showtime << std::fixed << std::setprecision(2) << std::scientific << SectionConsumedTime << "  "
                        << std::fixed << std::setprecision(2) << ((SectionConsumedTime/TotalConsumedTime)*100.0) << " %";
                Info::WriteStandard(Messages[it], showtime.str());
            }
        }
        Info::WriteLine("-");
        Info::WriteStandard("Total", to_string(TotalConsumedTime));
        Info::WriteLine("=");
    }
}