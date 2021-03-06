#ifndef INFO_H
#define INFO_H

//#include "Base/Includes.h"
#include <iomanip>
#include <sstream>
#include <string>
#include <type_traits>
#include <iostream>

namespace opensim
{
    class Settings;
    class dMatrix3x3;
    class dMatrix6x6;
    class dVector3;
    class dVector6;
    class Info
    {
        public:

        template <typename T>
        static inline std::string sgn(const T value)
        {
            return (value >= 0)?" ":"";
        }

        /// Alias for all "Standard" writing methods
        template <typename T>
        static void Write(const std::string Left, const  T Right,
                const bool showalways = true, const int precision = 2)
        {
            WriteStandard(Left, Right, showalways, precision);
        }
        /// Standard writing methods
        static void WriteStandard(const std::string Left, const std::string Right,
                const bool showalways = true, const int precision = 4);
        static void WriteStandard(const std::string Left, const dVector3 vec,
                const bool showalways = true, const int precision = 3);
        static void WriteStandard(const std::string Left, const dVector6 vec,
                const bool showalways = true, const int precision = 2);
        static void WriteStandard(const std::string Left, const dMatrix3x3 mat,
                const bool showalways = true, const int precision = 3);
        static void WriteStandard(const std::string Left, const dMatrix6x6 mat,
                const bool showalways = true, const int precision = 2);
        template <typename T>
        static void WriteStandard(const std::string Left, const T Right,
                const bool showalways = true, const int precision = 6)
        {
        #ifdef Silent
            if(showalways?true:false);
        #endif
            {
                if (std::is_floating_point<T>::value)
                {
                    std::cout << std::setfill(' ') << std::setw(ColumnLength)  << std::left << Left << ": ";
                    std::cout << std::setprecision(precision);
                    std::cout << std::scientific;
                    std::cout << sgn(Right) << Right << std::endl;
                }
                else if (std::is_signed<T>::value)
                {
                    std::cout << std::setfill(' ') << std::setw(ColumnLength)  << std::left << Left << ": ";
                    std::cout << sgn(Right) << Right << std::endl;
                }
                else
                {
                    std::cout << std::setfill(' ') << std::setw(ColumnLength)  << std::left << Left << ": ";
                    std::cout << " " << Right << std::endl;
                }
    
            }
        }

        static const size_t ColumnLength = 40;
        static const size_t LineLength   = 80;

        // Static console output methods
        static void WriteLine(const std::string LineType = "-",
                const bool showalways = true);
        static void WriteLineInsert(const std::string Insert,
                const std::string LineType = "-");
        static void WriteBlankLine(const bool showalways = true);
        static void WriteTimeStep(const int tStep, const int nSteps,
                const std::string Message = "");
        static void WriteTimeStep(const Settings& locSettings,
                const int tStep, const int nSteps,
                const std::string Message = "");
        static void WriteSimple(const std::string Message,
                const bool showalways = true);
        static void WriteStandardNarrow(const std::string Left,
                const std::string Right, const bool showalways = true);
        static void WriteCoordinate(const int x, const int y, const int z,
                const double dx, const bool showalways = true);
        static void WriteWarning(const std::string Message,
                const std::string Instance = "", const std::string Method = "",
                const bool showalways = true);
        static void WriteExit(const std::string Message,
                const std::string Instance = "", const std::string Method = "");
        static void WriteStartScreen(const bool showalways = true);
        static void PressEnterToContinue(void);
        static std::string GetStandard(const std::string Left,
                const std::string Right);
        static std::string GetStandardNarrow(const std::string Left,
                const std::string Right);

        template <typename T>
        static std::string to_string_with_precision(const T a_value,
                const int n = 6)
        {
            std::ostringstream out;
            out << std::setprecision(n) << a_value;
            return out.str();
        }
    };
}
#endif