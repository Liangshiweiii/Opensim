
#ifndef ELEMENT_H
#define ELEMENT_H

#include "Tools/Includes.h"

namespace opensim
{

struct Element                                                                  ///< Stores properties of chemcial element
{
    int             Index;                                                      ///< Element's index in components list.
    double          Value;                                                      ///< Element's amount
    double          Min;                                                        ///< This is the lower mole fraction limit of this component
    double          Max;                                                        ///< This is the highest mole fraction limit of this component
    std::string     Name;                                                       ///< Element's name (2 characters long max).
    bool            isVacancy;                                                  ///< Status of 'Component' or 'Vacancy'
    bool            isStoichiometric;                                           ///< Status if component has a solubility range in a specific phase (only to be used within phase storage)
    bool            Major;                                                      ///< Majority element in the sublattice
};

}
#endif//ELEMENT_H
