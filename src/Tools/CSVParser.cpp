

#include "Tools/CSVParser.h"
#include "Info.h"

namespace opensim
{
using namespace std;

void CSVParser::WriteHeader(const std::string fileName,
    const std::vector<std::string> headerArray, const std::string seperator)
{
    std::ofstream file(fileName.c_str());
    if (!file)
    {
        Info::WriteExit("File \"" + fileName + "\" could not be opened", "CSVParser", "WriteHeader()");
        exit(1);
    }
    if (headerArray.size() == 0)
    {
        Info::WriteWarning("Tried to write empty header array", "CSVParser", "WriteHeader()");
    }
    else
    {
        for (auto iter = headerArray.cbegin(); iter != headerArray.cend()-1; iter++)
        {
            file << *iter << seperator;
        }
        file << headerArray.back() << endl;
    }
    file.close();
}

void CSVParser::ClearContent(const std::string fileName)
{
    std::ofstream file(fileName.c_str(), ios_base::trunc);
    if (!file)
    {
        Info::WriteExit("File \"" + fileName + "\" could not be opened", "CSVParser", "ClearContent()");
        exit(1);
    }
    file.close();
}

void CSVParser::WriteData(const std::string fileName,
        const std::vector<int> dataArray, const std::string seperator)
{
    std::ofstream file(fileName.c_str(), ios_base::app);
    if (!file)
    {
        Info::WriteExit("File \"" + fileName + "\" could not be opened", "CSVParser", "WriteData()");
        exit(1);
    }
    if (dataArray.size() == 0)
    {
        Info::WriteWarning("Tried to write empty data array", "CSVParser", "WriteData()");
    }
    else
    {
        for (auto iter = dataArray.cbegin(); iter != dataArray.cend()-1; iter++)
        {
            file << *iter << seperator;
        }
        file << dataArray.back() << endl;
    }
    file.close();
}

void CSVParser::WriteData(const std::string fileName,
        const std::vector<double> dataArray, const std::string seperator)
{
    std::ofstream file(fileName.c_str(), ios_base::app);
    if (!file)
    {
        Info::WriteExit("File \"" + fileName + "\" could not be opened", "CSVParser", "WriteData()");
        exit(1);
    }
    if (dataArray.size() == 0)
    {
        Info::WriteWarning("Tried to write empty data array", "CSVParser", "WriteData()");
    }
    else
    {
        for (auto iter = dataArray.cbegin(); iter != dataArray.cend()-1; iter++)
        {
            file << *iter << seperator;
        }
        file << dataArray.back() << endl;
    }
    file.close();
}

void CSVParser::readFile(std::string& fileName,
    std::vector<std::vector<double>>& dataArray,
    std::vector<std::string>& headerArray, const std::string seperator)
{
    std::ifstream file(fileName.c_str());
    if (!file)
    {
        Info::WriteExit("File \"" + fileName + "\" could not be opened", "CSVParser", "readFile()");
        exit(1);
    };

    const char csep = seperator[0];
    bool headerFlag = false;
    std::string line;
    while (getline(file, line))
    {
        if (!headerFlag) //populate headers
        {
            std::stringstream lineStream(line);
            std::string cell;
            headerArray.clear();
            while (std::getline(lineStream, cell, csep))
            {
                headerArray.push_back(cell);
            }
            headerFlag = true;
        }
        else //populate data
        {
            std::vector<double> row_data;
            std::stringstream lineStream(line);
            std::string cell;
            row_data.clear();
            while (std::getline(lineStream, cell, csep))
            {
                row_data.push_back(atof(cell.c_str()));
            }
            dataArray.push_back(row_data);
        }
    }
    file.close();
}

void CSVParser::readFile(std::string& fileName,
    std::vector<std::vector<std::string>>& dataArray,
    std::vector<std::string>& headerArray, std::string seperator)
{
    std::ifstream file(fileName.c_str());
    if (!file)
    {
        Info::WriteExit("File \"" + fileName + "\" could not be opened", "CSVParser", "readFile()");
        exit(1);
    };

    const char csep = seperator[0];
    bool headerFlag = false;
    std::string line;
    while (getline(file, line))
    {
        if (!headerFlag) //populate headers
        {
            std::stringstream lineStream(line);
            std::string cell;
            headerArray.clear();
            while (std::getline(lineStream, cell, csep))
            {
                headerArray.push_back(cell);
            }
            headerFlag = true;
        }
        else //populate data
        {
            std::vector<std::string> row_data;
            std::stringstream lineStream(line);
            std::string cell;
            row_data.clear();
            while (std::getline(lineStream, cell, csep))
            {
                row_data.push_back(cell);
            }
            dataArray.push_back(row_data);
        }
    }
    file.close();
}

}// namespace opensim
