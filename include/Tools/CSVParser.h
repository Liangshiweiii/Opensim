#ifndef CSV_H
#define CSV_H

#include "Info.h"
#include "Tools/Includes.h"

namespace opensim
{

class CSVParser
{
 public:
    static void WriteHeader(const std::string fileName,
        const std::vector<std::string> headerArray, const std::string seperator = ",");
    static void ClearContent(const std::string fileName);
    static void WriteData(const std::string fileName,
        const std::vector<int> dataArray, const std::string seperator = ",");
    static void WriteData(const std::string fileName,
        const std::vector<double> dataArray, const std::string seperator = ",");
    static void readFile(std::string& fileName,
        std::vector<std::vector<double>>& dataArray,
        std::vector<std::string>& headerArray, const std::string seperator = ",");
    static void readFile(std::string& fileName,
        std::vector<std::vector<std::string>>& dataArray,
        std::vector<std::string>& headerArray, std::string seperator = ",");

 protected:
 private:
};

} // namespace opensim
#endif //CSV_H
