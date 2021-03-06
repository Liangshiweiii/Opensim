#include "Info.h"
#include "Tools/UserInterface.h"
#include "Settings.h"
#include "PhaseField.h"
#include "BoundaryConditions.h"
#include "Chemistry/ChemicalProperties.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "Mechanics/PlasticFlow/PlasticFlowCP.h"
#include "PhaseField.h"
#include "Settings.h"
#include "Temperatures.h"
#include "TextOutput.h"
#include <map>
#include "Compositions.h"

namespace opensim
{
using namespace std;

bool TextOutput::FileExists(string filename)
{
    /** This function checks the availability of the specified text file. If the
    file is nonexistent, locked or other problems prevent this function from
    opening the file, it will return the value 'false'.*/

    bool result = false;

    ifstream my_file(filename);

    if(my_file.good())
    result = true;

    return result;
}

void TextOutput::PhasePercent(ChemicalProperties& CP, PhaseField& Phi,
                                Settings& OPSettings, string filename,
                                double time)
{
    /** This function will create tabulated data on the volume percent of each
    thermodynamic phase defined in the chemical property input. Each time this
    function is called, a new row will be written in the specified file name
    for the current time step. If the file is not present it will be created,
    if the phase already exists, the new data will be appended!*/

    if (!FileExists(filename))
    {
        ofstream file(filename, ios::out);                                      //Create file and write header (or overwrite)
        file << "# time ";
        for(int idx = 0; idx < CP.Nphases; idx++)
        {
            file << CP.Phase[idx].Name << "  ";
        }
        file << endl;
        file.close();
    }

    ofstream file(filename, ios::app);                                          //Write content of the file
    vector<double> tPhaseFrac(CP.Nphases);
    double totalVolume = double(Phi.Nx * Phi.Ny * Phi.Nz);

    for(unsigned int idx = 0; idx < Phi.FieldsStatistics.size(); idx++)
    {
        tPhaseFrac[Phi.FieldsStatistics[idx].Phase]
                += Phi.FieldsStatistics[idx].Volume;
    }

    file << scientific << time << "   ";

    for(int idx = 0; idx < CP.Nphases; idx++)
    {
        file << scientific  <<  100.0*tPhaseFrac[idx]/totalVolume << "    ";
    }

    file << endl;
    file.close();                                                               //Close file
}

void TextOutput::AverageTemp(Temperature& Tx, string filename, double time)
{
    /** This function writes the average temperature of the system, each time it
    is called in a file.*/

    if (!FileExists(filename))
    {
        ofstream file(filename, ios::out);                                      //Create file and write header (or overwrite)
        file << "# Average Temperature in the system" << endl;
        file.close();
    }

    double average = 0.0;
    double points  = 0.0;
    STORAGE_LOOP_BEGIN(i,j,k,Tx.Tx,0)
    {
        points++;
        average += Tx.Tx(i,j,k);
    }
    STORAGE_LOOP_END

    if(points > 0.0)
    average /= points;

    ofstream file(filename, ios::app);                                          //Write content of the file

    file.precision(10);

    file << time << "    " << average << endl;

    file.close();
}

void TextOutput::GrainVolumes(PhaseField& Phi, string filename, double time)
{
    /** This function will write a list of all grain volumes in a file, each
    time this function is called in a new row. The first column is the time,
    the second row beginning with a "#" marks the thermodynamic index for each
    grain.*/

    int nPFs = Phi.FieldsStatistics.size();

    if (!FileExists(filename))
    {
        ofstream file(filename, ios::out);                                      //Create file and write header (or overwrite)
        file << "# Volume of each grain in the order of their field-index! "
             << "First row contains the simulation time. "
             << "The next line contains the thermodynamic phase field index "
             << "of each grain for identification" << endl;
        file << "#";
        for(int n = 0; n < nPFs; n++)
        {
            file << " " << Phi.FieldsStatistics[n].Phase;
        }
        file << endl;
        file.close();
    }

    map<int, double> AllGrains;

    STORAGE_LOOP_BEGIN(i,j,k,Phi.Fields,0)
    {
        for(auto n = Phi.Fields(i, j, k).cbegin();
                 n < Phi.Fields(i, j, k).cend(); ++n)
        {
            AllGrains[n->index] += n->value;
        }
    }
    STORAGE_LOOP_END

    ofstream file(filename, ios::app);                                          //Write content of the file

    file.precision(10);

    vector<double> AllGrainsSize(nPFs);
    AllGrainsSize.assign(nPFs, 0.0);
    for (map<int, double>::iterator it = AllGrains.begin();
            it != AllGrains.end(); ++it)
    {
        AllGrainsSize[it->first] = it->second;
    }

    file << time << " ";

    for (int i = 0; i < nPFs; i++)
    {
        file << AllGrainsSize[i] << " ";
    }
    file << endl;
    file.close();

    //This function was originally taken from PhaseField::WriteGrainsStatistics
}

void TextOutput::AverageStress(ElasticProperties& EP, string filename,
                                double timeOrStrain)
{
    /** This function will create tabulated data on the average stress in the
    system. Each time this function is called, a new row will be written in
    the specified file name for the current time step. If the file is not
    present it will be created, if the phase already exists, the new data
    will be appended!*/

    int Nx = EP.Nx;
    int Ny = EP.Ny;
    int Nz = EP.Nz;

    vStress avgStress;
    avgStress.set_to_zero();

    for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
    for (int k = 0; k < Nz; k++)
    {
        avgStress += EP.Stresses(i,j,k);
    }
    avgStress /= (Nx*Ny*Nz);

    if (!FileExists(filename))
    {
        ofstream file(filename, ios::out);
        file << "# time  "
             << "sigma_1  "
             << "sigma_2  "
             << "sigma_3  "
             << "sigma_4  "
             << "sigma_5  "
             << "sigma_6  "
             << "Pressure  "
             << "Mises  "
             << "sigma_norm";
        file << endl;
        file.close();
    }

    ofstream file(filename, ios::app);

    file << scientific << timeOrStrain << "   "
                << avgStress[0] << "  "
                << avgStress[1] << "  "
                << avgStress[2] << "  "
                << avgStress[3] << "  "
                << avgStress[4] << "  "
                << avgStress[5] << "  "
                << avgStress.Pressure() << "  "
                << avgStress.Mises() << "  "
                << avgStress.norm();
    file<< endl;
    file.close();
}

void TextOutput::AverageStrain(Storage3D<vStrain,0>& Strains, string filename,
                                double timeOrStrain)
{
    /** This function will create tabulated data on the average stress in the
    system. Each time this function is called, a new row will be written in
    the specified file name for the current time step. If the file is not
    present it will be created, if the phase already exists, the new data
    will be appended!*/

    int Nx = Strains.sizeX();
    int Ny = Strains.sizeY();
    int Nz = Strains.sizeZ();

    vStress avgStrain;
    avgStrain.set_to_zero();

    for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
    for (int k = 0; k < Nz; k++)
    {
        avgStrain += Strains(i,j,k);
    }
    avgStrain /=(Nx*Ny*Nz);

    if (!FileExists(filename))
    {
        ofstream file(filename, ios::out);
        file << "# time  "
             << "epsilon_1  "
             << "epsilon_2  "
             << "epsilon_3  "
             << "epsilon_4  "
             << "epsilon_5  "
             << "epsilon_6  ";
        file << endl;
        file.close();
    }

    ofstream file(filename, ios::app);

    file << scientific << timeOrStrain << "   "
                << avgStrain[0] << "  "
                << avgStrain[1] << "  "
                << avgStrain[2] << "  "
                << avgStrain[3] << "  "
                << avgStrain[4] << "  "
                << avgStrain[5];
    file<< endl;
    file.close();
}

void TextOutput::AverageDouble(Storage3D<double,0>& value, std::string filename,
                          double timeOrStrain)
{
    /** This function will create tabulated data on the average value in the
    system. Each time this function is called, a new row will be written in
    the specified file name for the current time step. If the file is not
    present it will be created, if the phase already exists, the new data
    will be appended!*/

    int Nx = value.sizeX();
    int Ny = value.sizeY();
    int Nz = value.sizeZ();

    double avgValue = 0;

    for (int i = 0; i < Nx; i++)
    for (int j = 0; j < Ny; j++)
    for (int k = 0; k < Nz; k++)
    {
        avgValue += value(i,j,k);
    }
    avgValue /=(Nx*Ny*Nz);

    if (!FileExists(filename))
    {
        ofstream file(filename, ios::out);
        file << "# time  "
             << "value";
        file << endl;
        file.close();
    }

    ofstream file(filename, ios::app);

    file << scientific << timeOrStrain << "   "
                << avgValue;
    file<< endl;
    file.close();
}

void TextOutput::AverageCRSS(PlasticFlowCP& PFCP, PhaseField& Phase, string filename,
                                double timeOrStrain)
{
    /** This function will create tabulated data on the CRSS in the
    system. Each time this function is called, a new row will be written in
    the specified file name for the current time step. If the file is not
    present it will be created, if the phase already exists, the new data
    will be appended!*/

    int Nx = PFCP.Nx;
    int Ny = PFCP.Ny;
    int Nz = PFCP.Nz;

    dVector<12> avgCRSS;
    avgCRSS.set_to_zero();
    double volume = Nx*Ny*Nz;
    STORAGE_LOOP_BEGIN(i,j,k,Phase.Fields,0)
    {
        if (Phase.Interface(i,j,k))
        {
            for(auto alpha = Phase.Fields(i, j, k).cbegin(); alpha < Phase.Fields(i, j, k).cend(); ++alpha)
            {
                int phaseIndex = alpha->index;
                for(int slipSys = 0; slipSys < 12; slipSys++)
                {
                    avgCRSS[slipSys] += PFCP.CRSS(i,j,k).get(phaseIndex, slipSys)*alpha->value/volume;
                }
            }
        }
        else
        {
            int phaseIndex = Phase.Fields(i,j,k).front().index;
            for(int slipSys = 0; slipSys < 12; slipSys++)
            {
                avgCRSS[slipSys] += PFCP.CRSS(i,j,k).get(phaseIndex, slipSys)/volume;
            }
        }
    }
    STORAGE_LOOP_END
    
    double avgtotal = 0.0;
    for(int slipSys = 0; slipSys < 12; slipSys++)
    {
        avgtotal += avgCRSS[slipSys]/12.0;
    }

    if (!FileExists(filename))
    {
        ofstream file(filename, ios::out);
        file << "# time  "
                << "total"
                << "ss_1"
                << "ss_2"
                << "ss_3"
                << "ss_4"
                << "ss_5"
                << "ss_6"
                << "ss_7"
                << "ss_8"
                << "ss_9"
                << "ss_10"
                << "ss_11"
                << "ss_12";
        file << endl;
        file.close();
    }

    ofstream file(filename, ios::app);

    file << scientific << timeOrStrain << "   "
            << avgtotal << "  "
            << avgCRSS[0] << "  "
            << avgCRSS[1] << "  "
            << avgCRSS[2] << "  "
            << avgCRSS[3] << "  "
            << avgCRSS[4] << "  "
            << avgCRSS[5] << "  "
            << avgCRSS[6] << "  "
            << avgCRSS[7] << "  "
            << avgCRSS[8] << "  "
            << avgCRSS[9] << "  "
            << avgCRSS[10] << "  "
            << avgCRSS[11];
    file<< endl;
    file.close();
}

void TextOutput::LineConcentration(Composition& Cx, ChemicalProperties& CP, PhaseField& Phi,
                                   string filename, double timestep,
                                   string type, string axis,
                                   int x, int y, int z)
{
    /** This function will create a separate file each time it is called, with
    a different file name according to the current time step. In this file, the
    tabulated composition data will be written down along a single line along a
    given axis with type = "X", "Y" or "Z". The line can be positioned with a
    point on it given with "x", "y" and "z". The output can be weight fraction
    with type = "WF", weight percent "WP", mole fraction "MF" and mole percent
    with "MP".*/

    int mode = -1;
    int dir = -1;
    int tablength = 16;
    std::stringstream ss;

    if(type == "WF")
    {
        mode = 0;
    }
    else if(type == "WP")
    {
        mode = 1;
    }
    else if(type == "MF")
    {
        mode = 2;
    }
    else if(type == "MP")
    {
        mode = 3;
    }
    else
    {
        cout << "Undefined type " << type
             << " in TextOutput::LineConcentration()! Chose \"WF\", \"WP\", "
             << "\"MF\" or \"MP\". Now exiting!" << endl;
        exit(1);
    }

    if(axis == "X")
    {
        dir = 0;
    }
    else if(axis == "Y")
    {
        dir = 1;
    }
    else if(axis == "Z")
    {
        dir = 2;
    }
    else
    {
        cout << "Undefined axis direction " << axis
             << " in TextOutput::LineConcentration()! Chose \"X\", \"Y\" or"
             << "\"Z\". Now exiting!" << endl;
        exit(1);
    }

    ss << fixed << std::setw(9) << std::setfill('0') << int(timestep);
    std::string s = ss.str();
    string name = filename + s + ".tab";

    ofstream file(name, ios::out);                                              //Create file and write header (or overwrite)

    file << std::setw(5) << std::setfill(' ') << "# x"
            << std::setw(5) << std::setfill(' ') << "y"
            << std::setw(5) << std::setfill(' ') << "z";

    for(int comp = 0; comp < CP.Ncomp; comp++)
    {
        string temp;
        string name(CP.Component[comp].Name);

        if(mode == 0)
        {
            temp = "wf. " + name;
        }
        else if(mode == 1)
        {
            temp = "wp. " + name;
        }
        else if(mode == 2)
        {
            temp = "mf. " + name;
        }
        else if(mode == 3)
        {
            temp = "mp. " + name;
        }

        file << std::setw(tablength) << std::setfill(' ') << temp;
    }

    file << endl;

    if(Cx.ConstituentFractions.size() > 0)
    CP.CalculateMoleFractions(Cx);

    if(dir == 0)                                                                //Write content of the file
    {
        for(int x = 0; x < CP.Nx; x++)
        {
            vector<double> Total(CP.Ncomp, 0.0);
            vector<double> WeightPercent(CP.Ncomp, 0.0);

            if(mode < 2)
            {
                WeightPercent = Cx.getWeightPercent(CP,Phi,x,y,z);
            }
            else
            {
                Total = CP.getTotalComposition(Cx,Phi,x,y,z);
            }

            file << fixed << std::setw(5) << std::setfill(' ') << x << ","
                          << std::setw(5) << std::setfill(' ') << y << ","
                          << std::setw(5) << std::setfill(' ') << z << ",";

            for(int comp = 0; comp < CP.Ncomp; comp++)
            {
                double val = 0.0;

                if(mode == 0)
                {//WF
                    val = WeightPercent[comp];
                }
                else if(mode == 1)
                {//WP
                    val = WeightPercent[comp]*100.;
                }
                else if(mode == 2)
                {//MF
                    val = Total[comp];
                }
                else if(mode == 3)
                {//MP
                    val = Total[comp]*100.;
                }

                file << scientific << std::setw(tablength) << std::setfill(' ')
                     << val << ",";
            }

            file << endl;
        }
    }
    else if(dir == 1)
    {
        for(int y = 0; y < CP.Ny; y++)
        {
            vector<double> Total(CP.Ncomp, 0.0);
            vector<double> WeightPercent(CP.Ncomp, 0.0);

            if(mode < 2)
            {
                WeightPercent = Cx.getWeightPercent(CP,Phi,x,y,z);
            }
            else
            {
                Total = CP.getTotalComposition(Cx,Phi,x,y,z);
            }

            file << fixed << std::setw(5) << std::setfill(' ') << x << ","
                          << std::setw(5) << std::setfill(' ') << y << ","
                          << std::setw(5) << std::setfill(' ') << z << ",";

            for(int comp = 0; comp < CP.Ncomp; comp++)
            {
                double val = 0.0;

                if(mode == 0)
                {//WF
                    val = WeightPercent[comp];
                }
                else if(mode == 1)
                {//WP
                    val = WeightPercent[comp]*100.;
                }
                else if(mode == 2)
                {//MF
                    val = Total[comp];
                }
                else if(mode == 3)
                {//MP
                    val = Total[comp]*100.;
                }

                file << scientific << std::setw(tablength) << std::setfill(' ')
                     << val << ",";
            }

            file << endl;
        }
    }
    else if(dir == 2)
    {
        for(int z = 0; z < CP.Nz; z++)
        {
            vector<double> Total(CP.Ncomp, 0.0);
            vector<double> WeightPercent(CP.Ncomp, 0.0);

            if(mode < 2)
            {
                WeightPercent = Cx.getWeightPercent(CP,Phi,x,y,z);
            }
            else
            {
                Total = CP.getTotalComposition(Cx,Phi,x,y,z);
            }

            file << fixed << std::setw(5) << std::setfill(' ') << x << ","
                          << std::setw(5) << std::setfill(' ') << y << ","
                          << std::setw(5) << std::setfill(' ') << z << ",";

            for(int comp = 0; comp < CP.Ncomp; comp++)
            {
                double val = 0.0;

                if(mode == 0)
                {//WF
                    val = WeightPercent[comp];
                }
                else if(mode == 1)
                {//WP
                    val = WeightPercent[comp]*100.;
                }
                else if(mode == 2)
                {//MF
                    val = Total[comp];
                }
                else if(mode == 3)
                {//MP
                    val = Total[comp]*100.;
                }

                file << scientific << std::setw(tablength) << std::setfill(' ')
                     << val << ",";
            }

            file << endl;
        }
    }

    file.close();                                                               //Close file
}

void TextOutput::LocalPhaseComposition(Composition& Cx, ChemicalProperties& CP,
                                       PhaseField& Phi, string filename,
                                       double time, int x, int y, int z)
{
    /** This function will create tabulated data on the phase composition of
    each thermodynamic phase defined in the chemical property input. Each time
    this function is called, a new row will be written in the specified file name
    for the current time step. If the file is not present it will be created,
    if the phase already exists, the new data will be appended!*/

    if (!FileExists(filename))
    {
        ofstream file(filename, ios::out);                                      //Create file and write header (or overwrite)
        file << "# time ";
        for(int alpha = 0; alpha < CP.Nphases; ++alpha)
        for(int comp = 0; comp < CP.Ncomp; comp++)
        {
            file << CP.Phase[alpha].Name << "_"
                 << CP.Phase[alpha].Component[comp].Name << "    ";
        }
        file << endl;
        file.close();
    }

    ofstream file(filename, ios::app);                                          //Write content of the file
    file << scientific << time << "   ";

    for(int alpha = 0; alpha < CP.Nphases; ++alpha)
    {
        if(Phi.Fractions(x,y,z)({alpha}) > 0.0)
        {
            Cx.Site2MoleFrac(CP.Phase[alpha],x,y,z);
            for(int comp = 0; comp < CP.Ncomp; comp++)
            {
                file << Cx.get_MF(alpha,comp,x,y,z) << "   ";
            }
        }
        else
        {            
            for(int comp = 0; comp < CP.Ncomp; comp++)
            file << 0.0 << "   ";
        }
    }

    file << endl;
    file.close();
}

}// namespace opensim
