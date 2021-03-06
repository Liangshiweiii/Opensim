#ifndef PHYSICALCONSTANTS_H
#define PHYSICALCONSTANTS_H

#include "Tools/Includes.h"

namespace opensim
{

class PhysicalConstants
{
 public:
    static double R;                                                            ///< Universal gas constant
    static constexpr double kBoltzmann = 1.3806488e-23;                         ///< Boltzmann constant [Joule/Kelvin]
    static constexpr double epsilon0   = 8.854187817620e-12;                    ///< Vacuum permittivity [ampere^2 second^4 Kg^(-1) meter^(-3)]
    static void ReadInput(std::string InputFileName);                           ///< Reads the physical constants values from the InputFile
};

}// namespace opensim
#endif
