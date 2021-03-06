#ifndef THERMODYNAMICPHASE_H
#define THERMODYNAMICPHASE_H

#include "Tools/Includes.h"
#include "Chemistry/SublatticeModel.h"

namespace opensim
{
class BoundaryConditions;
class ChemicalProperties;
class ElasticProperties;
class PhaseField;
class SublatticeModel;
class Settings;
class Composition;

struct PhaseInput
{
    int Number;
    std::vector<Element> Components;
    std::vector<std::string> RefElements;
};

/***/
class ThermodynamicPhase                                                        ///  Stores the properties of the thermodynamic phase
{
 public:
    //Constructors and structures
    ~ThermodynamicPhase(void);                                                  ///< Destructor
    ThermodynamicPhase(void);                                                   ///< Constructor
    ThermodynamicPhase(const ThermodynamicPhase& RHS)                           ///  Constructor
    {
        /**This constructor allows allocating and resizing of a vector of this
        phase, used in ChemicalProperties as std::vector<ThermodynamicPhase>*/
        Nx = RHS.Nx;            Ny = RHS.Ny;           Nz = RHS.Nz;
        Cmax = RHS.Cmax;      Cmin = RHS.Cmin;       Name = RHS.Name;
        Tmax = RHS.Tmax;      Tmin = RHS.Tmin;      Index = RHS.Index;
        Active = RHS.Active;                      Initial = RHS.Initial;
        ConsIndex = RHS.ConsIndex;             Sublattice = RHS.Sublattice;
        TotInitial = RHS.TotInitial;          MolarVolume = RHS.MolarVolume;
        Number = RHS.Number;
    }
    ThermodynamicPhase& operator=(const ThermodynamicPhase& RHS)                ///  Operator
    {
        /**This operator allows allocating and resizing of a vector of this
        phase, used in ChemicalProperties as std::vector<ThermodynamicPhase>*/
        if(this != &RHS)
        {
        Nx = RHS.Nx;            Ny = RHS.Ny;           Nz = RHS.Nz;
        Cmax = RHS.Cmax;      Cmin = RHS.Cmin;       Name = RHS.Name;
        Tmax = RHS.Tmax;      Tmin = RHS.Tmin;      Index = RHS.Index;
        Active = RHS.Active;                      Initial = RHS.Initial;
        ConsIndex = RHS.ConsIndex;             Sublattice = RHS.Sublattice;
        TotInitial = RHS.TotInitial;          MolarVolume = RHS.MolarVolume;
        Number = RHS.Number;
        }
        return *this;
    }
    struct CompositionLimits
    {
        double Min;
        double Max;
    };
    /*struct delta                                                                ///  Structure of in- and outgoing composition
    {
        //This structure stores all incoming and outgoing fluxes, used in the
        //diffusion and redistribution schemes
        double in;
        double out;
        bool inlimited;
        bool outlimited;
        double inlimit;
        double outlimit;
        delta& operator*= (double value)
        {
            in *= value;
            out *= value;
            inlimit *= value;
            outlimit *= value;
            return *this;
        }
        delta& operator+= (delta& locDelta)
        {
            in += locDelta.in;
            out += locDelta.out;
            inlimit += locDelta.in;
            outlimit += locDelta.out;
            return *this;
        }
        delta& operator-= (delta& locDelta)
        {
            in -= locDelta.in;
            out -= locDelta.out;
            inlimit -= locDelta.in;
            outlimit -= locDelta.out;
            return *this;
        }
    };*/

    double Nmoles(Composition& Cx, int x, int y, int z);                        ///< Returns the number of moles in the present phase
    //void SetBoundaryConditions(BoundaryConditions& BC);                       ///< Set boundary conditions for all storages
    void CalculateMoleFractions(Composition& Cx);                               ///< Calculate Molefractions for the whole domain
    void Initialize(Settings& locSettings, const int boundarysize, PhaseInput& phinput);///< Initialize this class
    int Nx, Ny, Nz;                                                             ///< Size of the inner calculation domain along X
    int Index;                                                                  ///< Phase index in e.g. TQ list of phases.
    int Number;                                                                 ///< Phase index as a position in the phase-vector
    int ReferenceElement;
    bool Active;                                                                ///< Indicates if current phase is considered in the current simulation
    double Tmin;                                                                ///< Indicates min temperature at which the phase is stable
    double Tmax;                                                                ///< Indicates max temperature at which the phase is stable
    Tensor<int,   1> ConsIndex;                                                 ///< Indicator where a sublattice starts and where it ends
    Tensor<double,1> MolarVolume;                                               ///< Molar volume for each phase //TODO: Storage3D
    Tensor<double,1> Cmin;                                                      ///< Minimum CompositionEXP value for each phase
    Tensor<double,1> Cmax;                                                      ///< Maximum CompositionEXP value for each phase
    Tensor<double,1> Initial;                                                   ///< Initial CompositionEXP of components in all phases
    Tensor<double,1> Initial2;                                                   ///< Initial CompositionEXP of components in all phases
    Tensor<double,1> Site2MoleFrac(Composition& Cx, int i, int j, int k);       ///< Will calculate Molefractions from Sitefractions in one point
    Tensor<double,1> InitialTracerPhaseComposition;                             ///< Tracer composition as given in input file
    std::string Name;                                                           ///< Phase name
    std::vector<double> TotInitial;                                             ///< Initial amount of each component
    std::vector<Element>         Component;                                     ///< List of all chemical components present in the phase
    std::vector<SublatticeModel> Sublattice;                                    ///< Sublattices composing the phase
    std::vector<double> yTotal(Composition& Cx, int x, int y, int z);           ///< Similar to NsitesWithI, but instead of a simple summation of all sites
    std::vector<double> unityWithoutVA(Composition& Cx, int x, int y, int z);   ///< Sum over all all site fractions in one sublattice (excluding VA)
    std::string thisclassname;                                                  ///< Name of this class
    //TODO: create constant values for and make them private:
    double Nsites;                                                              ///< Returns the sum of site coefficients of all sublattices
    int Ncons;                                                                  ///< Number of constituents in the present phase
    int Nsubs;                                                                  ///< Number of sublattices in the present phase
    int Ncomp;                                                                  ///< Number of components
    std::vector<int> SublIdx;                                                   ///< A vector which contains the first constituent in each sublattice
    std::vector<int> ConsIdx;                                                   ///< A vector that holds the Index of each Constituent
    std::vector<double> NsitesWithI;                                            ///< Summation of all sites of those sublattices, that include a certain element!
    std::vector<std::string> ConNames;                                          ///< A vector of strings with all the Names of each constituent
    //std::vector<Storage3D<delta,5> > SublDot;
    std::vector<CompositionLimits> SiteFractions;
    bool isInterstitial(int comp);                                              ///<
    bool isEndmember(int CIdxI);
    bool ispairwithVA(int CIdxI);
    bool isStoichiometric(int CIdxI);
    bool get_isStoichiometric(void);
    bool StoichiometricFlag;
    int Idx2Cons(int sub, int Idx);
    double Nm(Composition& Cx, int x, int y, int z);                            ///< obsolete
    std::vector<std::vector<double> > dMAdYi;
    bool AnalyticChemicalPotentials;
    /*SublatticeModel& Add_Sublattice(unsigned int n_sites)
    {
        SublatticeModel locSublattice(n_sites);
        unsigned int old_size = Constituent.size();
        Constituent.push_back(locElement);
        return Constituent[old_size];
    }*/

 protected:
 private:
 double getNsites(void);                                                        ///< Returns the sum of site coefficients of all sublattices
 int getNcons(void);                                                            ///< Number of constituents in the present phase
 int getNsubs(void);                                                            ///< Number of sublattices in the present phase
 int getNcomp(void);                                                            ///< Number of components
 std::vector<int> getSublIdx(void);                                             ///< A vector which contains the first constituent in each sublattice
 std::vector<int> getConsIdx(void);                                             ///< A vector that holds the Index of each Constituent
 std::vector<double> getNsitesWithI(void);                                      ///< Summation of all sites of those sublattices, that include a certain element!
 std::vector<std::string> getConNames(void);                                    ///< A vector of strings with all the Names of each constituent
 std::vector<std::vector<double> > getdMAdYi(void);
 bool getAnalyticChemicalPotentials(void);
};

}
#endif//THERMODYNAMICPHASE_H
