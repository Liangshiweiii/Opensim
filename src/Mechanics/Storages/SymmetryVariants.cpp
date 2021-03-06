/* This file is part of OPENPHASE BASE.
 * Copyright (c) 2020 OpenPhase Solutions GmbH, Bochum, Germany
 * For more details visit https://www.openphase-solutions.com
 * 
 * OPENPHASE BASE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *  
 * OPENPHASE BASE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *  
 * You should have received a copy of the GNU General Public License
 * along with OPENPHASE BASE.  If not, see <http://www.gnu.org/licenses/>. */

#include "Mechanics/Storages/SymmetryVariants.h"
#include "Settings.h"
#include "Info.h"
#include "Tools/UserInterface.h"

namespace opensim
{
using namespace std;

SymmetryVariants::SymmetryVariants(const Settings& locSettings, 
                                                 const std::string InputFileName)
{
    this->Initialize(locSettings);

    if(InputFileName == "default")
    {
        this->ReadInput();
    }
    else
    {
        this->ReadInput(InputFileName);
    }
}

void SymmetryVariants::Initialize(const Settings& locSettings)
{
    thisclassname = "SymmetryVariants";
    //DefaultInputFileName = ProjectInputDir + "SymmetryVariantsInput.opi";

    Nphases = locSettings.Nphases;
    TransformationMatrices.resize(Nphases);
    set = false;
    Info::WriteStandard("SymmetryVariants", "Initialized");
}

void SymmetryVariants::ReadInput(string InputFileName)
{
    Info::WriteLine();
    Info::WriteLineInsert("Crystal symmetry transformation matrices");
    Info::WriteStandard("Source", InputFileName);

    fstream inp(InputFileName.c_str(), ios::in);

    if (!inp)
    {
        Info::WriteExit("File " + InputFileName + " could not be opened", thisclassname, "ReadInput()");
        exit(1);
    };
    
    int moduleLocation = UserInterface::FindModuleLocation(inp, thisclassname);
    
    // Reading number of variants for each phase
    for(unsigned int pIndex = 0; pIndex < Nphases; pIndex++)
    {
        stringstream converter;
        converter << "Nvariants" << "_" << pIndex;
        int locNvariants = UserInterface::ReadParameterI(inp, moduleLocation, converter.str(),true,0);
        TransformationMatrices[pIndex].resize(locNvariants);
        for (int n = 0; n < locNvariants; n++)
        {
            TransformationMatrices[pIndex][n].set_to_unity();            
        }
    }
    // Reading variants for each phase
    for(unsigned int pIndex = 0; pIndex < Nphases; pIndex++)
    for(unsigned int n = 0; n < TransformationMatrices[pIndex].size(); n++)
    {
        stringstream converter;
        converter << "V" << "_" << pIndex << "_" << n;
        UserInterface::FindParameterLocation(inp, moduleLocation, converter.str());
        TransformationMatrices[pIndex][n].read(inp);
    }   
    set = true;
    inp.close();

    stringstream outstream;

    Info::WriteBlankLine();
    for (unsigned int pIndex = 0; pIndex < Nphases; pIndex++)
    for (unsigned int vIndex = 0; vIndex < TransformationMatrices[pIndex].size(); vIndex++)
    {
        outstream << "T(" << pIndex << ", " << vIndex << "):\n" << TransformationMatrices[pIndex][vIndex].print() << endl;
    }
    Info::WriteSimple(outstream.str());
    Info::WriteLine();
}

}// namespace openphase
