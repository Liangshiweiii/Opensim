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

#ifndef VTK_H
#define VTK_H

#include "Tools/Includes.h"
namespace opensim
{

const int PDScalars = 1; // Plain scalar entries
const int PDVectors = 2; // 3 component vector
const int PDNormals = 3; // 3 component vector
const int PDTensors = 4; // 3x3 matrix

class VTK                                                                      /// static class to write VTK data in xml format
{
 public:
    // Partial VTK output methods
    static void WriteHeader(std::stringstream& buffer,
                                         const int Nx,
                                         const int Ny,
                                         const int Nz);
    static void WriteBeginPointData(std::stringstream& buffer,
                const std::vector<int> PointDataTypes);
    static void WriteEndPointData(std::stringstream& buffer);

    static void WriteCoordinates(std::stringstream& buffer,
                                         const int Nx,
                                         const int Ny,
                                         const int Nz);
    static void WriteCoordinatesParallel(std::stringstream& buffer,
                                         const int Nx, const int offx,
                                         const int Ny, const int offy,
                                         const int Nz, const int offz);
    static void WriteToFile(std::stringstream& buffer,
                                  const std::string FileName,
                                  const int tStep);

    // Partial VTK output methods for scalar values
    static void WriteScalar(std::stringstream& buffer,
                                const Storage3D<int, 0>& OutputData,
                                  const std::string VariableName);
    static void WriteScalar(std::stringstream& buffer,
                                const Storage3D<double, 0>& OutputData,
                                  const std::string VariableName);
    // Partial VTK output methods for vector containg scalars (arbitrary components)
    static void WriteVectorComponent(std::stringstream& buffer,
                                const Storage3D<int, 1>& OutputData,
                                  const std::string VariableName,
                                  const int component);
    static void WriteVectorComponent(std::stringstream& buffer,
                                const Storage3D<double, 1>& OutputData,
                                  const std::string VariableName,
                                  const int component);
    // Partial VTK output methods for OpenPhase dVector3 (three components)
    static void WriteVector(std::stringstream& buffer,
                                  const Storage3D<dVector3,0>& OutputData,
                                  const std::string VariableName);
    static void WriteVector(std::stringstream& buffer,
                                  const Storage3D<dVector3,1>& OutputData,
                                  const std::string VariableName,
                                  const int component);
    static void WriteVectorAsVector(std::stringstream& buffer,
                                  const Storage3D<dVector3,0>& OutputData,
                                  const std::string VariableName);
    static void WriteVectorAsVector(std::stringstream& buffer,
                                  const Storage3D<dVector3,1>& OutputData,
                                  const std::string VariableName,
                                  const int component);
    // Partial VTK output methods for OpenPhase dMatrix3x3
    static void WriteMatrix(std::stringstream& buffer,
                                  const Storage3D<dMatrix3x3,0>& OutputData,
                                  const std::string VariableName);
    // Partial VTK output methods for symmetric tensors (six components)
    static void WriteSymmetricTensor(std::stringstream& buffer,
                                  const Storage3D<dVector6,0>& OutputData,
                                  const std::string VariableName);
    static void WriteSymmetricTensor(std::stringstream& buffer,
                                  const Storage3D<vStrain,0>& OutputData,
                                  const std::string VariableName);
    static void WriteSymmetricTensor(std::stringstream& buffer,
                                  const Storage3D<vStress,0>& OutputData,
                                  const std::string VariableName);
    // Partial VTK output methods for tensors (arbitrary components)
    static void WriteTensor(std::stringstream& buffer,
                                  const Storage3D<int,2>& OutputData,
                                  const std::string VariableName,
                                  const int componentA, const int componentB);
    static void WriteTensor(std::stringstream& buffer,
                                  const Storage3D<double,2>& OutputData,
                                  const std::string VariableName,
                                  const int componentA, const int componentB);
    // Partial VTK output methods for rotations
    static void WriteQuaternion(std::stringstream& buffer,
                                  const Storage3D<Quaternion,0>& OutputData,
                                  const std::string VariableName);
    static void WriteAngles(std::stringstream& buffer,
                                const Storage3D<EulerAngles,0>& OutputData,
                                const std::string VariableName,
                                EulerConvention EConvention);

};
}// namespace opensim
#endif
