#ifndef SERIAL
    #include <omp.h>
#else // overloading several OpenMP methods used in OpenPhase
    inline int omp_get_num_threads()
    {
        return 1;
    }
    inline int omp_get_max_threads()
    {
        return 1;
    }
    inline int omp_get_thread_num()
    {
        return 0;
    }
    inline void omp_set_num_threads(int){};
#endif

#include <vector>
#include <array>
#include <algorithm>
#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>
#include <cmath>
#include <complex>
#include <cfloat>
#include <limits>
#include <errno.h>
#include <signal.h>
#include <stdexcept>
#include <random>
#include <ciso646>

#ifndef INCLUDES_H
#define INCLUDES_H

namespace opensim
{
inline void simulation_end()
{
#ifdef _WIN32
    std::getchar();                                                             ///< Prevents closing of the terminal window at the end of the simulation
#endif // _WIN32
}

#ifdef _WIN32
    const std::string dirSeparator = "\\";                                      ///< Windows style directory separator
#else
    const std::string dirSeparator = "/";                                       ///< Unix/Linux style directory separator
#endif

//#ifndef Pi
const double Pi = 3.14159265358979323846;                                       ///< Pi constant value
//#endif

const std::complex< double > I(0.0, 1.0);                                       ///< sqrt(-1) declaration


typedef double (*ThermoDynFnctPtr)(std::vector<double>& nu, double T);          ///< Multi component thermodynamic function, takes a vector of density/concentration values as input
typedef double (*CPAnalyticFnctPtr)(std::vector<double>& nu, double T, double P);///< Multi component thermodynamic function, takes a vector of density/concentration values as input

// Special keywords definition
const int Default = 0;
const int Random = 1;

//Aggregate state markers
const int Solid  = 0;
const int Liquid = 1;
const int Gas    = 3;

const int SolidSolid   = 0; /*Solid + Solid*/
const int LiquidLiquid = 2; /*Liquid + Liquid*/
const int GasGas       = 6; /*Gas + Gas*/
const int SolidLiquid  = 1; /*Solid + Liquid*/
const int SolidGas     = 3; /*Solid + Gas*/
const int LiquidGas    = 4; /*Liquid + Gas*/

// Schemes for advection
const int Upwind      = 0;
const int LaxWendroff = 1;
const int VanLeer     = 2;
const int Superbee    = 3;

const double LaplacianStencil27[3][3][3] = {{{1.0/30.0,   1.0/10.0, 1.0/30.0},
                                             {1.0/10.0,   7.0/15.0, 1.0/10.0},
                                             {1.0/30.0,   1.0/10.0, 1.0/30.0}},

                                            {{1.0/10.0,   7.0/15.0, 1.0/10.0},
                                             {7.0/15.0, -64.0/15.0, 7.0/15.0},
                                             {1.0/10.0,   7.0/15.0, 1.0/10.0}},

                                            {{1.0/30.0,   1.0/10.0, 1.0/30.0},
                                             {1.0/10.0,   7.0/15.0, 1.0/10.0},
                                             {1.0/30.0,   1.0/10.0, 1.0/30.0}}};///< 27 point Laplacian stencil by Spotz and Carey (1995)

const double LaplacianStencil27b[3][3][3] ={{{1.0/48.0,   1.0/8.0, 1.0/48.0},
                                             {1.0/8.0,   5.0/12.0, 1.0/8.0},
                                             {1.0/48.0,   1.0/8.0, 1.0/48.0}},

                                            {{1.0/8.0,   5.0/12.0, 1.0/8.0},
                                             {5.0/12.0, -25.0/6.0, 5.0/12.0},
                                             {1.0/8.0,   5.0/12.0, 1.0/8.0}},

                                            {{1.0/48.0,   1.0/8.0, 1.0/48.0},
                                             {1.0/8.0,   5.0/12.0, 1.0/8.0},
                                             {1.0/48.0,   1.0/8.0, 1.0/48.0}}}; ///< 27 point Laplacian stencil by Dave Hale

const double LaplacianStencil19[3][3][3] = {{{1.0/30.0,   1.0/10.0, 1.0/30.0},
                                             {1.0/10.0,   7.0/15.0, 1.0/10.0},
                                             {1.0/30.0,   1.0/10.0, 1.0/30.0}},

                                            {{1.0/10.0,   7.0/15.0, 1.0/10.0},
                                             {7.0/15.0, -64.0/15.0, 7.0/15.0},
                                             {1.0/10.0,   7.0/15.0, 1.0/10.0}},

                                            {{1.0/30.0,   1.0/10.0, 1.0/30.0},
                                             {1.0/10.0,   7.0/15.0, 1.0/10.0},
                                             {1.0/30.0,   1.0/10.0, 1.0/30.0}}};///< 19 point Laplacian stencil (not yet implemented)

const double LaplacianStencil7[3][3][3] =  {{{0.0,   0.0, 0.0},
                                             {0.0,   1.0, 0.0},
                                             {0.0,   0.0, 0.0}},

                                            {{0.0,   1.0, 0.0},
                                             {1.0,  -6.0, 1.0},
                                             {0.0,   1.0, 0.0}},

                                            {{0.0,   0.0, 0.0},
                                             {0.0,   1.0, 0.0},
                                             {0.0,   0.0, 0.0}}};               ///< 7 point Laplacian stencil

template<typename T>
inline void ignore_result(T /* unused result */) {}                             ///< Allows to suppress warnings on unused function return values

const std::string DefaultInputFileName = "ProjectInput.opi";
const std::string ProjectInputDir = "ProjectInput" + dirSeparator;
const std::string VTKDir = "VTK" + dirSeparator;
const std::string RawDataDir = "RawData" + dirSeparator;

}// namespace opensim
#endif

#include "Tools/Macros.h"
#include "Tools/Storages.h"
#include "Tools/OPObjects.h"
#include "PhysicalConstants.h"