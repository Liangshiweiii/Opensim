#include "Tools/UserInterface.h"
#include "Info.h"
namespace opensim
{
using namespace std;

// Auxiliary functions
string removeEntireWhiteSpaces(const string InString)
{
    string StringTemp = InString;
    StringTemp.erase(std::remove_if( StringTemp.begin(), StringTemp.end(),
     [](char c){ return (c =='\r' || c =='\t' || c == ' ' || c == '\n');}), StringTemp.end() );
     return StringTemp;
}

string removeLeadingTrailingWhiteSpaces(const string InString)
{
    std::string whitespaces (" \t");
    size_t first = InString.find_first_not_of(whitespaces);
    size_t last = InString.find_last_not_of(whitespaces);
    if (first == last) return "";
    return InString.substr(first, (last-first+1));
}

int StringToInt(const std::string& str, const std::string name)
{
    int ivar = 0;

    try
    {
        ivar = std::stoi(str);
    }
    catch (const std::invalid_argument&)
    {
        Info::WriteExit("Argument for $" + name + " is invalid", "UserInterface", "ReadParameterI()");
        //throw;
        exit(3);
    }
    return ivar;
}

double StringToDouble(const std::string& str, const std::string name)
{
    double dvar = 0.0;

    try
    {
        dvar = std::stod(str);
    }
    catch (const std::invalid_argument&)
    {
        Info::WriteExit("Argument for $" + name + " is invalid", "UserInterface", "ReadParameterD()");
        //throw;
        exit(3);
    }
    return dvar;
}
// Auxiliary functions end

string UserInterface::MakeFileName(string Directory, string NameBase, int Index, string FileExtension)
{
    stringstream converter;
    converter << Index;
    string Count = converter.str();
    int n = Count.size();
    string Number;
    for(int i = 0; i < 8 - n; i++) Number.append("0");
    Number.append(Count);

    return Directory + NameBase + Number + FileExtension;
}

int UserInterface::FindModuleLocation(fstream& Inp, const string module)
{
    // Returns the location of the first line in the module specified in the parameters
    Inp.seekg(0, ios::beg);
    while (Inp.good())
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '@');
        if (Inp.eof()) break;
        Inp >> ReadKey;

        if (!ReadKey.compare(module))
        {
            streampos temp = Inp.tellg();
            //temp += 1;
            return temp;
        }
    }
    Inp.seekg(0, ios::beg);
    return 0;
}

int UserInterface::FindParameter(fstream& Inp,
    const int location, const string Key)
{
    streampos curPos(location+1);
    Inp.seekg(curPos);

    bool found = false;
    int keyLocation = 0;
    while (!Inp.eof() and !found)
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            keyLocation = -1;
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                keyLocation = -1;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            keyLocation = Inp.tellg();
            keyLocation -= (ReadKey.size() + 1);
            found = true;
            break;
        }
    }
    Inp.clear();
    return keyLocation;
}

int UserInterface::FindParameterLocation(fstream& Inp,
    const int location, const string Key)
{
    streampos curPos(location+1);
    Inp.seekg(curPos);

    bool found = false;
    int keyLocation = 0;
    while (!Inp.eof() and !found)
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            keyLocation = -1;
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                keyLocation = -1;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            //getline(Inp, tmp, ':');
            keyLocation = Inp.tellg();
            //keyLocation -= (ReadKey.size() + 1);
            found = true;
            break;
        }
    }
    Inp.clear();
    return keyLocation;
}

double UserInterface::ReadParameterD(fstream& Inp, int currentLocation, string Key,
    const bool mandatory, const double defaultval)
{
    // Note the following optional arguments:
    // verbose: if true, found variable is printed to terminal, standard = true
    // mandatory: if parameter not found, program will stop (standard = 1 = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    double ReturnValue = defaultval;

    while (!Inp.eof())
    {
        string tmp;
        string tmp2;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            getline(Inp, tmp2);
            tmp = removeLeadingTrailingWhiteSpaces(tmp);
            tmp2 = removeEntireWhiteSpaces(tmp2);

            // Transfer string tmp2 to double (and check if this is possible at all)
            ReturnValue = StringToDouble(tmp2, tmp);

            // Checks if comment is given, otherwise use "Key" as output
            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::Write(Key, ReturnValue, true, 5);
            }
            else
            {
                Info::Write(tmp, ReturnValue, true, 5);
            }
            found = true;
            break;
        }
    }

    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead
            ReturnValue = defaultval;
        }
    }
    Inp.clear();
    return ReturnValue;
}

string UserInterface::ReadParameterS(fstream& Inp, int currentLocation, string Key,
    const bool mandatory, const string defaultval)
{
    // Note the following optional arguments:
    // verbose: if true, found variable is printed to terminal, standard = true
    // mandatory: if parameter not found, program will stop (standard = 1 = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    string ReturnValue;
    while (!Inp.eof())
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            // Yet problem with blanks in name
            getline(Inp, tmp, ':');
            Inp >> ReturnValue;
            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::WriteStandard(Key, ReturnValue);
            }
            else
            {
                tmp.erase(0, tmp.find_first_not_of(" \t"));
                tmp.erase(tmp.find_last_not_of(" \t") + 1, tmp.size());
                Info::WriteStandard(tmp, ReturnValue);
            }
            found = true;
            break;
        }
    }

    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead
            ReturnValue = defaultval;
        }
    }

    Inp.clear();
    return ReturnValue;
}

int UserInterface::ReadParameterB(fstream& Inp, int currentLocation, const string Key,
    const bool mandatory, const string defaultval)
{
    // Note the following optional arguments:
    // verbose: if true, found variable is printed to terminal, standard = true
    // mandatory: if parameter not found, program will stop (standard = 1 = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    int ReturnValue = 0;
    string Answer;
    while (!Inp.eof())
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            getline(Inp, Answer);
            Answer.erase(0, Answer.find_first_not_of("YesNo") + 1);
            if (Answer.find_last_of("YesNo") + 1 != Answer.size())
            {
                Answer.erase(Answer.find_last_of("YesNo") + 1, Answer.size());
            }
            //Inp >> ReturnValue;
            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::WriteStandard(Key, Answer);
            }
            else
            {
                tmp.erase(0, tmp.find_first_not_of(" \t"));
                tmp.erase(tmp.find_last_not_of(" \t") + 1, tmp.size());
                Info::WriteStandard(tmp, Answer);
            }
            found = true;
            break;
        }
    }
    string tmp;
    for (unsigned int i = 0; i < Answer.size(); ++i)
        if (Answer[i] != ' ') tmp.push_back(Answer[i]);

    if (found and !tmp.compare("Yes")) ReturnValue = 1;
    if (found and !tmp.compare("No")) ReturnValue = 0;
    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            ReturnValue = 0;
        }
    }
    Inp.clear();
    return ReturnValue;
}

string UserInterface::ReadParameterF(fstream& Inp, int currentLocation, const string Key,
    const bool mandatory, const string defaultval)
{
    // Note the following optional arguments:
    // verbose: if true, found variable is printed to terminal, standard = true
    // mandatory: if parameter not found, program will stop (standard = 1 = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval

    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    string ReturnValue;
    string tFileName;
    while (!Inp.eof())
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            getline(Inp, tFileName);
            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::WriteStandard(Key, tFileName);
            }
            else
            {
                tmp.erase(0, tmp.find_first_not_of(" \t"));
                tFileName.erase(0, tFileName.find_first_not_of(" \t"));
                Info::WriteStandard(tmp, tFileName);
            }
            found = true;
            break;
        }
    }
    string tmp;
    for (unsigned int i = 0; i < tFileName.size(); i++)
        if (tFileName[i] != ' ') ReturnValue.push_back(tFileName[i]);

    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead
            ReturnValue = defaultval;
        }
    }

    Inp.clear();
    return ReturnValue;
}

int UserInterface::ReadParameterI(fstream& Inp, int currentLocation, string Key,
    const bool mandatory, int const defaultval)
{
    // Note the following optional arguments:
    // verbose: if true, found variable is printed to terminal, standard = true
    // mandatory: if parameter not found, program will stop (standard = true = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    // currentLocation: current location reading is done at, if new section is found/file ends,
    // -1 is returned; if not, updated location is returned
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    int ReturnValue = defaultval;
    while (!Inp.eof())
    {
        string tmp;
        string tmp2;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if (!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            getline(Inp, tmp2);
            tmp = removeLeadingTrailingWhiteSpaces(tmp);
            tmp2 = removeEntireWhiteSpaces(tmp2);

            // Transfer string tmp2 to int (and check if this is possible at all)
            ReturnValue = StringToInt(tmp2, tmp);

            if (tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::Write(Key, ReturnValue);
            }
            else
            {
                Info::Write(tmp, ReturnValue);
            }
            found = true;
            break;
        }
    }

    if (not found)
    {
        if (mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead
            ReturnValue = defaultval;
        }
    }
    Inp.clear();
    return ReturnValue;
}

char UserInterface::ReadParameterC(fstream& Inp, int currentLocation, const string Key,
                              const bool mandatory, const char defaultval)
{
    // Note the following optional arguments:
    // verbose: if true, found variable is printed to terminal, standard = true
    // mandatory: if parameter not found, program will stop (standard = 1 = mandatory)
    // defaultval: if parameter not found and not mandatory, method returns defaultval
    streampos curPos(currentLocation);
    Inp.seekg(curPos);
    bool found = false;
    char ReturnValue;
    while (!Inp.eof())
    {
        string tmp;
        string ReadKey;

        getline(Inp, tmp, '$');
        if (Inp.eof())
        {
            break;
        }

        bool endReading = false;
        for (unsigned int i = 0; i < tmp.size(); i++)
        {
            if (tmp.at(i) == '@')
            {
                endReading = true;
                break;
            }
        }

        if (endReading)
            break;

        Inp >> ReadKey;

        if(!ReadKey.compare(Key))
        {
            getline(Inp, tmp, ':');
            Inp >> ReturnValue;
            if(tmp.find_first_not_of(' ') == std::string::npos)
            {
                Info::Write(Key, tmp);
            }
            else
            {
                tmp.erase(0,tmp.find_first_not_of(" \t"));
                Info::Write(tmp, ReturnValue);
            }
            found = true;
            break;
        }
    }

    if(not found)
    {
        if(mandatory)
        {
            Info::WriteExit("Parameter for Key \"" + Key + "\" does not exist in the input file.");
            exit(1);
        }
        else
        {
            // Return defaultval instead ( = 0 if not defined explicitly)
            ReturnValue = defaultval;
        }
    }
    Inp.clear();
    return ReturnValue;
}

}//namespace opensim
