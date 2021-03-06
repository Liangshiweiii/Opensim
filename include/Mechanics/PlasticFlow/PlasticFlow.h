

#ifndef PLASTICFLOW_H
#define PLASTICFLOW_H

#include "Tools/Includes.h"

namespace opensim
{
class Settings;
class PhaseField;
class ElasticProperties;
class Composition;
class Temperature;
class Settings;
class Orientations;
class BoundaryConditions;

class PlasticFlow : public OPObject
{
public:
    virtual ~PlasticFlow(){};

    virtual void Initialize(Settings& locSettings) = 0;
    using OPObject::ReadInput;
    virtual void ReadInput(std::string InputFileName) = 0;

    virtual void SetBoundaryConditions(const BoundaryConditions& BC) = 0;
    virtual void Remesh(int newNx, int newNy, int newNz, BoundaryConditions& BC)
                                                                            = 0;

    virtual double Solve(std::vector<OPObject*> Objects, bool verbose) = 0;

    double dt;

protected:
private:

};
}// namespace opensim

#endif
