#ifndef TIMEINFO_H
#define TIMEINFO_H

#include "Tools/Includes.h"
#include "Settings.h"
#include "Info.h"
namespace opensim
{

class Settings;
class TimeInfo
{
 public:

    TimeInfo(){};
    TimeInfo(const Settings& locSettings, const std::string Name = "NONE");
    ~TimeInfo(void);
    void Initialize(const Settings& locSettings, const std::string Name);
    void SetStart(void);
    void SetTimeStamp(const std::string Message);
    void SkipToHere(void);
    void Reset(void);
    void PrintWallClockSummary(void) const;
    void PrintCPUClockSummary(void) const
    {
        Info::WriteWarning("TimeInfo::PrintCPUClockSummary() not yet implemented.",
                "TimeInfo", "PrintCPUClockSummary()");
    };
    void PrintFullSummary(void) const
    {
        Info::WriteWarning("TimeInfo::PrintFullSummary() not yet implemented.",
                "TimeInfo", "PrintFullSummary()");
    };


 protected:
 private:
    std::string thisclassname;
    std::vector<clock_t> TimeCollector;
    std::vector<std::string> Messages;

    std::string TimerName;
};
}// namespace opensim
#endif