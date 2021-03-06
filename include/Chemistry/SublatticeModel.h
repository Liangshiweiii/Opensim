
#ifndef SUBLATTICEMODEL_H
#define SUBLATTICEMODEL_H

#include "Tools/Includes.h"
#include "Chemistry/Element.h"

namespace opensim
{

class SublatticeModel                                                           ///< Stores the sublattice properties
{
 public:
    int Ncons;
    bool hasVacancies;
    void Initialize(void);
    SublatticeModel(double n_sites, unsigned int n_cons = 0);                   ///< Constructor
    //bool hasVacancies();                                                      ///< Returns boolean values, if sublattice has vacancies
    bool isElementPresent(int i);                                               ///< Returns boolean values, if constituent is present on sublattice
    //int Ncons();                                                              ///< Returns the number of constituents on the sublattice
    std::vector<Element> Constituent;                                           ///< Constituents composing the sublattice
    double Site;                                                                ///< Sublattice contribution to the chemical formular unit
    int Reference;

    /*Element& Add_Constituent(unsigned int index, double value, char* name, bool isVA)
    {
        Element locElement(index, value, name, isVA);
        unsigned int old_size = Constituent.size();
        Constituent.push_back(locElement);    

        return Constituent[old_size];
    }*/

 protected:
 private:
    int getNcons();
    bool gethasVacancies();
};

}
#endif//SUBLATTICEMODEL_H
