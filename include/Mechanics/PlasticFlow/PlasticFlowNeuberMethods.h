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

#ifndef PLASTICFLOWNEUBERMETHODS_H
#define PLASTICFLOWNEUBERMETHODS_H

#include "Tools/Includes.h"

namespace opensim
{
	class ElasticProperties;

	class PlasticFlowNeuberMethods
	{
	public:
		static void getNeuberData(vStress StressIn, vStrain StrainIn, vStress& StressOut, vStrain& StrainOut);
		static void getNeuberDataRO(vStress StressIn, vStrain StrainIn, vStress& StressOut, vStrain& StrainOut);
		static vStrain getNeuberStrains(vStrain StrainIn);
		static vStress getNeuberStresses(vStress StressIn);

		static void writeNeuberStressVTK(ElasticProperties& EP, const int tStep);
		static void writeNeuberStrainVTK(ElasticProperties& EP, const int tStep);

	protected:
		static void WriteStressesVTKData(ElasticProperties& EP, std::stringstream& buffer);
		static void WriteStrainsVTKData(ElasticProperties& EP, std::stringstream& buffer);

	private:

	};
}// namespace openphase

#endif

