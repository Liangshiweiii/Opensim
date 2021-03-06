#ifndef CRYSTALLOGRAPHY_H
#define CRYSTALLOGRAPHY_H

#include "Tools/Includes.h"

namespace opensim
{

class Crystallography : OPObject                                                /// Stores the mass densities of phases in each point.
{
 public:

    Crystallography(){};
    Crystallography(const Settings& locSettings, std::string InputFileName = "NONE");
    void Initialize(const Settings& locSettings);                               /// Initializes storages, sets internal variables.

    Storage<dMatrix3x3> SymmetriesCubic;                                        /// Stores cubic symmetry operations

 protected:
 private:
};

} // namespace openphase

#endif