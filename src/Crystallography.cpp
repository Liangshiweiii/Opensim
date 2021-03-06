#include "Crystallography.h"
#include "Info.h"

namespace opensim
{
using namespace std;

Crystallography::Crystallography(const Settings& locSettings, const std::string InputFileName)
{
    this->Initialize(locSettings);

//    if(InputFileName == "NONE")
//    {
//        this->ReadInput();
//    }
//    else
//    {
//        this->ReadInput(InputFileName);
//    }
}

void Crystallography::Initialize(const Settings& locSettings)
{
    thisclassname = "Crystallography";

    // Symmetry matrices for cubic lattice. Taken from Kocks, Tome and Wenk
    // 'Texture and Anisotropy', Cambridge University Press 1998

    dMatrix3x3 SCubic1  { 1, 0, 0,   0, 1, 0,   0, 0, 1 };
    dMatrix3x3 SCubic2  { 0, 0,-1,   0,-1, 0,  -1, 0, 0 };
    dMatrix3x3 SCubic3  { 0, 0, 1,   1, 0, 0,   0, 1, 0 };
    dMatrix3x3 SCubic4  { 0, 0, 1,   0,-1, 0,   1, 0, 0 };
    dMatrix3x3 SCubic5  { 0, 1, 0,   0, 0, 1,   1, 0, 0 };
    dMatrix3x3 SCubic6  { 0, 0, 1,   0, 1, 0,  -1, 0, 0 };

    dMatrix3x3 SCubic7  { 0,-1, 0,   0, 0, 1,  -1, 0, 0 };
    dMatrix3x3 SCubic8  { 0, 0,-1,   0, 1, 0,   1, 0, 0 };
    dMatrix3x3 SCubic9  { 0,-1, 0,   0, 0,-1,   1, 0, 0 };
    dMatrix3x3 SCubic10 {-1, 0, 0,   0, 0,-1,   0,-1, 0 };
    dMatrix3x3 SCubic11 { 0, 1, 0,   0, 0,-1,  -1, 0, 0 };
    dMatrix3x3 SCubic12 { 1, 0, 0,   0, 0,-1,   0, 1, 0 };

    dMatrix3x3 SCubic13 { 0, 0,-1,   1, 0, 0,   0,-1, 0 };
    dMatrix3x3 SCubic14 { 1, 0, 0,   0, 0, 1,   0,-1, 0 };
    dMatrix3x3 SCubic15 { 0, 0,-1,  -1, 0, 0,   0, 1, 0 };
    dMatrix3x3 SCubic16 {-1, 0, 0,   0, 0, 1,   0, 1, 0 };
    dMatrix3x3 SCubic17 { 0, 0, 1,  -1, 0, 0,   0,-1, 0 };
    dMatrix3x3 SCubic18 { 0,-1, 0,  -1, 0, 0,   0, 0,-1 };

    dMatrix3x3 SCubic19 {-1, 0, 0,   0, 1, 0,   0, 0,-1 };
    dMatrix3x3 SCubic20 { 0, 1, 0,  -1, 0, 0,   0, 0, 1 };
    dMatrix3x3 SCubic21 {-1, 0, 0,   0,-1, 0,   0, 0, 1 };
    dMatrix3x3 SCubic22 { 0, 1, 0,   1, 0, 0,   0, 0,-1 };
    dMatrix3x3 SCubic23 { 1, 0, 0,   0,-1, 0,   0, 0,-1 };
    dMatrix3x3 SCubic24 { 0,-1, 0,   1, 0, 0,   0, 0, 1 };

    SymmetriesCubic.Allocate(24);

    SymmetriesCubic[0] = SCubic1;
    SymmetriesCubic[1] = SCubic2;
    SymmetriesCubic[2] = SCubic3;
    SymmetriesCubic[3] = SCubic4;
    SymmetriesCubic[4] = SCubic5;
    SymmetriesCubic[5] = SCubic6;
    SymmetriesCubic[6] = SCubic7;
    SymmetriesCubic[7] = SCubic8;
    SymmetriesCubic[8] = SCubic9;
    SymmetriesCubic[9] = SCubic10;
    SymmetriesCubic[10] = SCubic11;
    SymmetriesCubic[11] = SCubic12;
    SymmetriesCubic[12] = SCubic13;
    SymmetriesCubic[13] = SCubic14;
    SymmetriesCubic[14] = SCubic15;
    SymmetriesCubic[15] = SCubic16;
    SymmetriesCubic[16] = SCubic17;
    SymmetriesCubic[17] = SCubic18;
    SymmetriesCubic[18] = SCubic19;
    SymmetriesCubic[19] = SCubic20;
    SymmetriesCubic[20] = SCubic21;
    SymmetriesCubic[21] = SCubic22;
    SymmetriesCubic[22] = SCubic23;
    SymmetriesCubic[23] = SCubic24;

    initialized = true;
    Info::WriteStandard(thisclassname, "Initialized");
}

}// namespace opensim