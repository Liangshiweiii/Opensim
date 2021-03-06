#include "Info.h"
#include "Settings.h"

namespace opensim
{
    using namespace std;

    void Info::WriteTimeStep(const int tStep, const int nSteps, const string Message)
    {
        time_t rawtime;
        time(&rawtime);

        string sTime;
        string thetime(ctime( &rawtime ));
        sTime = thetime.substr(0, thetime.length() -1);

        WriteLine("_");
        cout << setfill(' ') << setw(ColumnLength)  << left <<
                "TimeStep" << " " << to_string(tStep) + "/" + to_string(nSteps) << endl;
        cout << setfill(' ') << setw(ColumnLength)  << left <<
                "Time" << " " << sTime << endl;
        if (Message != "")
        {
            cout << Message << endl;
        }
        WriteLine("_");
    }

    void Info::WriteTimeStep(const Settings& locSettings, const int tStep,
            const int nSteps, const string Message)
    {
        if (!(tStep%locSettings.tScreenWrite))
        {
            WriteTimeStep(tStep, nSteps, Message);
        }
    }

    void Info::WriteLine(const string LineType, const bool showalways)
    {
    #ifndef SILENT
        //char cstr[LineType.length()+1];
        //strcpy(cstr, LineType.c_str());
        cout << setfill(LineType[0]) << setw(LineLength) << "" << endl;
    #else
        if(showalways)
        {
            //char cstr[LineType.length()+1];
            //strcpy(cstr, LineType.c_str());
            cout << setfill(LineType[0]) << setw(LineLength) << "" << endl;
        }
    #endif
    }

    void Info::WriteLineInsert(const string Insert, const string LineType)
    {
        //char cstr[LineType.length()+1];
        //strcpy(cstr, LineType.c_str());
        cout << LineType[0] << LineType[0] << "< " << Insert << " >" << setfill(LineType[0]) <<
                setw(std::max<int>(0,LineLength-Insert.size()-6)) << "" << endl;
    }

    void Info::WriteBlankLine(const bool showalways)
    {
    #ifndef SILENT
        cout << endl;
    #else
        if(showalways)
        {
            cout << endl;
        }
    #endif
    }

    void Info::WriteSimple(const string Message, const bool showalways)
    {
    #ifndef SILENT
        cout << Message << endl;
    #else
        if(showalways)
        {
            cout << Message << endl;
        }
    #endif
    }

    void Info::WriteStandard(const string Left, const string Right,
            const bool showalways, const int precision)
    {
    #ifndef SILENT
        cout << /*setfill(' ') << setw(ColumnLength)  <<*/ left << Left << ":  " << Right << endl;
    #else
        if(showalways)
        {
            cout << /*setfill(' ') << setw(ColumnLength)  <<*/ left << Left << ":  " << Right << endl;
        }
    #endif
    }

    void Info::WriteStandard(const std::string Left, const dVector3 vec,
            const bool showalways, const int precision)
    {
    #ifdef Silent
        if(showalways?true:false);
    #endif
        {
            std::cout << /*setfill(' ') << setw(ColumnLength) <<*/ left << Left << ": ";
            std::cout << setprecision(precision);
            std::cout << std::scientific;

            std::cout << "|"
                << sgn(vec[0]) << vec[0] << " "
                << sgn(vec[1]) << vec[1] << " "
                << sgn(vec[2]) << vec[2] << " |\n";
        }
    }

    void Info::WriteStandard(const std::string Left, const dVector6 vec,
            const bool showalways, const int precision)
    {
    #ifdef Silent
        if(showalways?true:false);
    #endif
        {
            std::cout << setfill(' ') << setw(ColumnLength)  << left << Left << ": ";

            std::cout << setprecision(precision);
            std::cout << std::scientific;

            std::cout << "|"
                << sgn(vec[0]) << vec[0] << " "
                << sgn(vec[1]) << vec[1] << " "
                << sgn(vec[2]) << vec[2] << " "
                << sgn(vec[3]) << vec[3] << " "
                << sgn(vec[4]) << vec[4] << " "
                << sgn(vec[5]) << vec[5] << " |\n";
        }
    }

    void Info::WriteStandard(const std::string Left, const dMatrix3x3 mat,
            const bool showalways, const int precision)
    {
    #ifdef Silent
        if(showalways?true:false);
    #endif
        {
            std::cout << setfill(' ') << setw(ColumnLength)  << left << Left << ": ";
            for (unsigned int i = 0;  i < 3; ++i)
            {
                std::cout << setprecision(precision);
                std::cout << std::scientific;

                std::cout << "|"
                    << sgn(mat(i,0)) << mat(i,0) << " "
                    << sgn(mat(i,1)) << mat(i,1) << " "
                    << sgn(mat(i,2)) << mat(i,2) << "  |\n";
            }
        }
    }

    void Info::WriteStandard(const std::string Left, const dMatrix6x6 mat,
            const bool showalways, const int precision)
    {
    #ifdef Silent
        if(showalways?true:false);
    #endif
        {
            std::cout << setfill(' ') << setw(ColumnLength)  << left << Left << ": ";
            for (unsigned int i = 0;  i < 6; ++i)
            {
                std::cout << setprecision(precision);
                std::cout << std::scientific;

                std::cout << "|"
                    << sgn(mat(i,0)) << mat(i,0) << " "
                    << sgn(mat(i,1)) << mat(i,1) << " "
                    << sgn(mat(i,2)) << mat(i,2) << " "
                    << sgn(mat(i,3)) << mat(i,3) << " "
                    << sgn(mat(i,4)) << mat(i,4) << " "
                    << sgn(mat(i,5)) << mat(i,5) << " |\n";
            }
        }
    }

    void Info::WriteStandardNarrow(const string Left, const string Right, const bool showalways)
    {
    #ifndef SILENT
        cout << setfill(' ') << setw(20)  << left << Left << ": " << Right << endl;
    #else
        if(showalways)
        {
            cout << setfill(' ') << setw(20)  << left << Left << ": " << Right << endl;
        }
    #endif
    }

    string Info::GetStandard(const string Left, const string Right)
    {
        stringstream ss;
        ss << setfill(' ') << setw(ColumnLength)  << left << Left << ": " << setprecision(8) << Right << endl;
        string returnStandard = ss.str();

        return returnStandard;
    }

    string Info::GetStandardNarrow(const string Left, const string Right)
    {
        stringstream ss;
        ss << setfill(' ') << setw(20)  << left << Left << ": " << setprecision(8) << Right << endl;
        string returnStandard = ss.str();

        return returnStandard;
    }

    void Info::WriteCoordinate(const int x, const int y, const int z, const double dx,
                                            const bool showalways)
    {
    #ifndef SILENT
        cout << "Point:      (" << x << ", " << y << ", " << z << ")" << endl;
        cout << "Coordinate: (" << x*dx << ", " << y*dx << ", " << z*dx << ")" << endl;
    #else
        if(showalways)
        {
            cout << "Point:      (" << x << ", " << y << ", " << z << ")" << endl;
            cout << "Coordinate: (" << x*dx << ", " << y*dx << ", " << z*dx << ")" << endl;
        }
    #endif
    }
    void Info::WriteWarning(const string Message, const string Instance,
                            const string Method, const bool showalways)
    {
    #ifndef SILENT
        string thisInstance = Instance;
        WriteLine("~");
        if (Method != "")
        {
            thisInstance += "::" + Method;
        }
        cout << setfill(' ') << setw(10)  << left << "Warning: "  << thisInstance << endl
                                                << "          " << Message  << endl;
        WriteLine("~");
    #else
        if(showalways)
        {
            string thisInstance = Instance;
            WriteLine("~");
            if (Method != "")
            {
                thisInstance += "::" + Method;
            }
            cout << setfill(' ') << setw(10)  << left << "Warning: "  << thisInstance << endl
                                                    << "          " << Message  << endl;
            WriteLine("~");
        }
    #endif
    }

    void Info::WriteExit(const string Message, const string Instance, const string Method)
    {
        string thisInstance = Instance;
        time_t rawtime;
        time(&rawtime);

        string sTime;
        string thetime(ctime( &rawtime ));
        sTime = thetime.substr(0, thetime.length() -1);

        cerr << endl;
        WriteLine("*");
        if (Method != "")
        {
            thisInstance += "::" + Method;
        }
        cerr << setfill(' ') << setw(10)  << left << "Calculation terminated!" << endl
                            << setw(10)  << left << "Instance:" << thisInstance << endl
                            << setw(10)  << left << "Reason: " << Message << endl
                            << setw(10)  << left << "Time: " << sTime << endl;
        WriteLine("*");
        cerr << endl;
    }

    void Info::WriteStartScreen(const bool showalways)
    {
    #ifndef SILENT
        cout << endl;
        WriteLine(">");
        cout << "  OpenPhase" << endl << endl
            << "  An open-source phase-field simulation library" << endl
            << "  www.openphase.de" << endl << endl
            << "  Core development: " << endl
            << "  Interdisciplinary Centre for Advanced Materials Simulation (ICAMS) " << endl
            << "  Ruhr University Bochum, Germany " << endl
            << "  2009 - 2018" << endl << endl
            << "  Licensed under GNU GPLv3" << endl;
        WriteLine("<");
        cout << endl;
    #else
        if(showalways)
        {
            cout << endl;
            WriteLine(">");
            cout << "  OpenPhase" << endl << endl
                << "  An open-source phase-field simulation library" << endl
                << "  www.openphase.de" << endl << endl
                << "  Core development: " << endl
                << "  Interdisciplinary Centre for Advanced Materials Simulation (ICAMS) " << endl
                << "  Ruhr University Bochum, Germany " << endl
                << "  2009 - 2018" << endl << endl
                << "  Licensed under GNU GPLv3" << endl;
            WriteLine("<");
            cout << endl;
        }
    #endif
    }

    void Info::PressEnterToContinue(void)
    {
        cout << "Press ENTER to continue... " << flush;
        cin.ignore( numeric_limits <streamsize> ::max(), '\n' );
    }
}