#ifndef USERINTERFACE_H
#define USERINTERFACE_H

#include "Tools/Includes.h"

namespace opensim
{
class UserInterface
{
 public:
    static std::string MakeFileName(std::string Directory,
                                      std::string NameBase,
                                      int Index,
                                      std::string FileExtension);               /// Creates full filename string using the location directory, name base, running index and file extension

    // Methods to read input parameters from OpenPhase input files.
    // Searches for $KEY in the entire file using the following syntax:
    // $KEY    commment    :   value
    
    static int FindParameter(std::fstream& Inp,                                 /// Checks if Key is present in the specified module and returns location right before the key
		                     const int location, std::string Key);
		                         
	static int FindParameterLocation(std::fstream& Inp,                         /// Checks if Key is present in the specified module and returns location right after the key
		                     const int location, std::string Key);
		                         
	static int FindModuleLocation(std::fstream& Inp,							/// Returns location of the module
								const std::string module);

	static double ReadParameterD(std::fstream& Inp,							    /// Read double precision floating point parameter value
		int currentLocation,
		const std::string Key,
		const bool mandatory = true,
		const double defaultval = 0.0);

	static int ReadParameterI(std::fstream& Inp,							    /// Read integer parameter value
		int currentLocation,
		std::string Key,
		const bool mandatory = true,
		int const defaultval = 0);

	static std::string ReadParameterS(std::fstream& Inp,                        /// Read next string parameter value
		int currentLocation,
		const std::string Key,
		const bool mandatory = true,
		const std::string defaultval = "NN");

	static int ReadParameterB(std::fstream& Inp,                                /// Read boolean parameter value
		int currentLocation,
		std::string Key,
		const bool mandatory = true,
		const std::string defaultval = "No");

	static std::string ReadParameterF(std::fstream& Inp,					    /// Read filename string (assumes no spaces in the filename)
		int currentLocation,
		const std::string Key, 
		const bool mandatory = true, 
		const std::string defaultval = "NN");

	static char ReadParameterC(std::fstream& Inp,                               /// Read char parameter value
		int currentLocation,
		const std::string Key,
		const bool mandatory = true,
		const char defaultval = 'X');
};
}// namespace opensim
#endif
