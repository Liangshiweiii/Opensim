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

#include "Mechanics/PlasticFlow/PlasticFlowNeuberMethods.h"
#include "Mechanics/Storages/ElasticProperties.h"
#include "VTK.h"

namespace opensim
{
	using namespace std;

	void PlasticFlowNeuberMethods::getNeuberData(vStress StressIn, vStrain StrainIn, vStress& StressOut, vStrain& StrainOut)
	{
	    for(int n = 0; n < 6; n++)
	    {
	        double strain_sign = 1;
	        if(StrainIn[n] < 0.0) strain_sign = -1;
	        double stress_sign = 1;
	        if(StressIn[n] < 0.0) stress_sign = -1;
	        double energy = fabs(StressIn[n] * StrainIn[n]*1.0e-6);	    
	        double newStrain = pow(energy/1217.0, 2.0/3.0);
	        double newStress = 1217.0*sqrt(newStrain);
	        StressOut[n] = newStress*stress_sign*1.0e6;
	        StrainOut[n] = newStrain*strain_sign;
	    }
	}
    void PlasticFlowNeuberMethods::getNeuberDataRO(vStress StressIn, vStrain StrainIn, vStress& StressOut, vStrain& StrainOut)
	{
	    double Modulus = 160.0e9; // Young's modulus
	    double sigma0 = 200.0e6;  // Yield stress
	    double alpha = 3.0; // Materials parameter (exponemnt n = 3)
	    for(int n = 0; n < 6; n++)
	    {
	        double strain_sign = 1;
	        if(StrainIn[n] < 0.0) strain_sign = -1;
	        double stress_sign = 1;
	        if(StressIn[n] < 0.0) stress_sign = -1;
	        double energy = fabs(StressIn[n] * StrainIn[n]);	    

	        double newStress = sqrt(-0.5*sigma0*sigma0/alpha + 0.5*sqrt((4.0*Modulus*energy*alpha*sigma0*sigma0 + pow(sigma0,4))/pow(alpha,2)));
	        double newStrain = newStress/Modulus*(1.0 + alpha*pow(newStress/sigma0,2));
	        StressOut[n] = newStress*stress_sign;
	        StrainOut[n] = newStrain*strain_sign;
	    }
	}
	vStrain PlasticFlowNeuberMethods::getNeuberStrains(vStrain StrainIn)
	{
		return StrainIn;			//Put Neuber equations here
	}

	vStress PlasticFlowNeuberMethods::getNeuberStresses(vStress StressIn)
	{
		return StressIn;			//Put Neuber equations here
	}

	void PlasticFlowNeuberMethods::writeNeuberStressVTK(ElasticProperties& EP, const int tStep)
	{
		stringstream buffer;
		vector<int> DataTypes{PDScalars};

		VTK::WriteHeader(buffer, EP.Nx, EP.Ny, EP.Nz);
		VTK::WriteBeginPointData(buffer, DataTypes);
		{
			PlasticFlowNeuberMethods::WriteStressesVTKData(EP, buffer);
		}
		VTK::WriteEndPointData(buffer);
		VTK::WriteCoordinates(buffer, EP.Nx, EP.Ny, EP.Nz);
		VTK::WriteToFile(buffer, "NeuberStresses", tStep);
	}

	void PlasticFlowNeuberMethods::writeNeuberStrainVTK(ElasticProperties& EP, const int tStep)
	{
		stringstream buffer;
		vector<int> DataTypes{ PDScalars };
		VTK::WriteHeader(buffer, EP.Nx, EP.Ny, EP.Nz);
		VTK::WriteBeginPointData(buffer, DataTypes);
		{
			PlasticFlowNeuberMethods::WriteStrainsVTKData(EP, buffer);
		}
		VTK::WriteEndPointData(buffer);
		VTK::WriteCoordinates(buffer, EP.Nx, EP.Ny, EP.Nz);
		VTK::WriteToFile(buffer, "NeuberStrains", tStep);
	}

	void PlasticFlowNeuberMethods::WriteStressesVTKData(ElasticProperties& EP, stringstream& buffer)
	{
		vector<long int> dimV{EP.Stresses.sizeX(), EP.Stresses.sizeY(), EP.Stresses.sizeZ()};
		vector<int> compV{ 0, 1, 2, 3, 4, 5 };
		vector<string> compNameV{ "xx", "yy", "zz", "yz", "xz", "xy" };

		for (auto it = compV.cbegin(); it != compV.cend(); ++it)
		{
			string compname = "\"NeuberStress_" + compNameV[*it] + "\" ";
			buffer << "<DataArray type = \"Float64\" Name = " << compname <<
				"NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < EP.Stresses.sizeZ(); ++k)
				for (int j = 0; j < EP.Stresses.sizeY(); ++j)
					for (int i = 0; i < EP.Stresses.sizeX(); ++i)
					{
						buffer << getNeuberStresses(EP.Stresses(i, j, k))[*it] << "\n";
					}
			buffer << "</DataArray>" << endl;
		}

		buffer << "<DataArray type = \"Float64\" Name = \"" << "Pressure" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < EP.Stresses.sizeZ(); ++k)
			for (int j = 0; j < EP.Stresses.sizeY(); ++j)
				for (int i = 0; i < EP.Stresses.sizeX(); ++i)
				{
					buffer << getNeuberStresses(EP.Stresses(i, j, k)).Pressure() << endl;
				}
		buffer << "</DataArray>" << endl;

		buffer << "<DataArray type = \"Float64\" Name = \"" << "vMises" <<
			"\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
		for (int k = 0; k < EP.Stresses.sizeZ(); ++k)
			for (int j = 0; j < EP.Stresses.sizeY(); ++j)
				for (int i = 0; i < EP.Stresses.sizeX(); ++i)
				{
					buffer << getNeuberStresses(EP.Stresses(i, j, k)).Mises() << endl;
				}
		buffer << "</DataArray>" << endl;
	}

	void PlasticFlowNeuberMethods::WriteStrainsVTKData(ElasticProperties& EP, stringstream& buffer)
	{
		vector<int> compV{ 0, 1, 2, 3, 4, 5 };
		vector<string> compNameV{ "1", "2", "3", "4", "5", "6" };

		for (auto it = compV.cbegin(); it != compV.cend(); ++it)
		{
			buffer << "<DataArray type = \"Float64\" Name = \""
				<< "E_" << compNameV[*it]
				<< "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
			for (int k = 0; k < EP.Nz; ++k)
				for (int j = 0; j < EP.Ny; ++j)
					for (int i = 0; i < EP.Nx; ++i)
					{
						buffer << getNeuberStrains(EP.Strains(i, j, k))[*it] << endl;
					}
			buffer << "</DataArray>" << endl;
		}
	}
}
