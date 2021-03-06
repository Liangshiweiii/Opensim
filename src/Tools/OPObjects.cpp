#include "Tools/OPObjects.h"
#include "Info.h"
#include "Settings.h"

namespace opensim
{
bool OPObject::IsInitialized(void) const
{
    return initialized;
}

bool OPObject::IsNotInitialized(void) const
{
    return !initialized;
}

OPObject::~OPObject(void)
{
    if(this->thisclassname.size() !=0)
    {
        Info::WriteStandard(this->thisclassname, "Exited normally");
    }
    else
    {
        Info::WriteWarning("Uninitialized object exited normally",
                           "OPObject","~OPObject");
    }
}

void OPObject::ReadInput(void)
{
    if(this->thisclassname != "")
    {
        this->ReadInput(DefaultInputFileName);
    }
    else
    {
        Info::WriteExit("Attempt to use an uninitialized object!", "OPObject",
                                                                   "ReadInput");
        exit(13);
    }
}

}// namespace opensim