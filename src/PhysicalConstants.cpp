#include "Info.h"
#include "PhysicalConstants.h"
#include "Tools/UserInterface.h"

namespace opensim
{

using namespace std;
/*************************************************************************/
double PhysicalConstants::R = 8.314;                                            // in Jouls/(mole*K)

void PhysicalConstants::ReadInput(string InputFileName)
{
    fstream inp(InputFileName.c_str(), ios::in);

    if (!inp)
    {
        Info::WriteExit("File " + InputFileName + " could not be opened","PhysicalConstants");
        exit(1);
    };

    int moduleLocation = UserInterface::FindModuleLocation(inp, "PhysicalConstants");
    R = UserInterface::ReadParameterD(inp,moduleLocation, string("R"));

    inp.close();
}

}// namespace opensim
