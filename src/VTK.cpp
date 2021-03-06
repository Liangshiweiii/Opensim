#include "Settings.h"
#include "VTK.h"
#include "Tools/UserInterface.h"
#include "Info.h"
namespace opensim
{

using namespace std;

vector<string> voigtcon {"xx", "yy", "zz", "yz", "xz", "xy"};
vector<string> matrixcon  {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"};

void VTK::WriteHeader(stringstream& buffer, const int Nx, const int Ny, const int Nz)
{
    buffer << "<?xml version= \"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>" << endl;
    buffer << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
    buffer << "<StructuredGrid WholeExtent=\""
             << 0 << " " << Nx-1 << " "
             << 0 << " " << Ny-1 << " "
             << 0 << " " << Nz-1 << "\"> " << endl;
}

void VTK::WriteBeginPointData(stringstream& buffer, const vector<int> PointDataTypes)
{
    // Method requires vector that contains the types of the given point data.
    // Each type has to be given once. The following numbers (or corresponding constants)
    // can be used:
    // PDScalars = 1
    // PDVectors = 2
    // PDNormals = 3
    // PDTensors = 4

    // Example vector initialization for an output with scalar and tensorial data:
    // ' std::vector<int> DataTypes {PDScalars, PDTensors}; '

    buffer << "<PointData ";

    for(auto it = PointDataTypes.cbegin(); it != PointDataTypes.cend(); ++it)
    {
        switch(*it)
        {
            case PDScalars:
            {
                buffer << " Scalars= \"ScalarData\"";
                break;
            }
            case PDVectors:
            {
                buffer << " Vectors= \"VectorData\"";
                break;
            }
            case PDNormals:
            {
                buffer << " Normals= \"NormalData\"";
                break;
            }
            case PDTensors:
            {
                buffer << " Tensors= \"TensorData\"";
                break;
            }
            default:
            {
                break;
            }
        }
    }
    buffer << ">" << endl;
}

void VTK::WriteEndPointData(stringstream& buffer)
{
    buffer << "</PointData>" << endl;
}

void VTK::WriteCoordinates(stringstream& buffer, const int Nx, const int Ny, const int Nz)
{
    buffer << "<Points>" << endl;
    buffer << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        buffer << i << " " << j << " " << k << "\n";
    }
    buffer << "</DataArray>" << endl;
    buffer << "</Points>" << endl;
}

void VTK::WriteCoordinatesParallel(stringstream& buffer, const int Nx,
        const int offx, const int Ny, const int offy, const int Nz, const int offz)
{
    buffer << "<Points>" << endl;
    buffer << "<DataArray type = \"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    for(int k = 0; k < Nz; ++k)
    for(int j = 0; j < Ny; ++j)
    for(int i = 0; i < Nx; ++i)
    {
        buffer << i+offx << " " << j+offy << " " << k+offz << "\n";
    }
    buffer << "</DataArray>" << endl;
    buffer << "</Points>" << endl;
}

void VTK::WriteToFile(stringstream& buffer,
        const std::string FileName, const int tStep)
{
    buffer << "</StructuredGrid>" << endl;
    buffer << "</VTKFile>" << endl;
    string completeFileName = UserInterface::MakeFileName(VTKDir,
                                             FileName + "_", tStep, ".vts");

    ofstream vtk_file(completeFileName.c_str());
    vtk_file << buffer.rdbuf();
    vtk_file.close();
}

// Write Point Data methods

void VTK::WriteScalar(stringstream& buffer, const Storage3D<int, 0>& OutputData,
                         const string VariableName)
{
    // Corresponding PointDataTypes: PDScalars = 1

    buffer << "<DataArray type = \"Int64\" Name = \"" << VariableName <<
                "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for(int k = 0; k < OutputData.sizeZ(); ++k)
    for(int j = 0; j < OutputData.sizeY(); ++j)
    for(int i = 0; i < OutputData.sizeX(); ++i)
    {
        buffer << OutputData(i,j,k) << "\n";
    }
    buffer << "</DataArray>" << endl;
}

void VTK::WriteScalar(stringstream& buffer, const Storage3D<double, 0>& OutputData,
                         const string VariableName)
{
    // Corresponding PointDataTypes: PDScalars = 1

    buffer << "<DataArray type = \"Float64\" Name = \"" << VariableName <<
                "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for(int k = 0; k < OutputData.sizeZ(); ++k)
    for(int j = 0; j < OutputData.sizeY(); ++j)
    for(int i = 0; i < OutputData.sizeX(); ++i)
    {
        buffer << OutputData(i,j,k) << "\n";
    }
    buffer << "</DataArray>" << endl;
}

void VTK::WriteVectorComponent(stringstream& buffer, const Storage3D<int, 1>& OutputData,
                         const string VariableName, const int component)
{
    // Corresponding PointDataTypes: PDScalars = 1

    buffer << "<DataArray type = \"Int64\" Name = \"" << VariableName
           << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for(int k = 0; k < OutputData.sizeZ(); ++k)
    for(int j = 0; j < OutputData.sizeY(); ++j)
    for(int i = 0; i < OutputData.sizeX(); ++i)
    {
        buffer << OutputData(i,j,k)({component}) << "\n";
    }
    buffer << "</DataArray>" << endl;
}

void VTK::WriteVectorComponent(stringstream& buffer, const Storage3D<double, 1>& OutputData,
                         const string VariableName, const int component)
{
    // Corresponding PointDataTypes: PDScalars = 1

    buffer << "<DataArray type = \"Float64\" Name = \""
            << VariableName
            << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for(int k = 0; k < OutputData.sizeZ(); ++k)
    for(int j = 0; j < OutputData.sizeY(); ++j)
    for(int i = 0; i < OutputData.sizeX(); ++i)
    {
        buffer << OutputData(i,j,k)({component}) << "\n";
    }
    buffer << "</DataArray>" << endl;
}

void VTK::WriteVector(stringstream& buffer,
        const Storage3D<dVector3,0>& OutputData, const string VariableName)
{
    // Corresponding PointDataTypes: PDScalars = 1

    for (int n = 0; n < 3; n++)
    {
        buffer << "<DataArray type = \"Float64\" Name = \""
               << VariableName + "_" << to_string(n)
               << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < OutputData.sizeZ(); ++k)
        for(int j = 0; j < OutputData.sizeY(); ++j)
        for(int i = 0; i < OutputData.sizeX(); ++i)
        {
            buffer << OutputData(i,j,k)[n] << "\n";
        }
        buffer << "</DataArray>" << endl;
    }
}

void VTK::WriteVector(stringstream& buffer,
        const Storage3D<dVector3,1>& OutputData, const string VariableName,
        const int component)
{
    // Corresponding PointDataTypes: PDScalars = 1

    for (int n  = 0; n  < 3; n++)
    {
        buffer << "<DataArray type = \"Float64\" Name = \""
               << VariableName << "_" << n
               << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < OutputData.sizeZ(); ++k)
        for(int j = 0; j < OutputData.sizeY(); ++j)
        for(int i = 0; i < OutputData.sizeX(); ++i)
        {
            buffer << OutputData(i,j,k)({component})[n] << "\n";
        }
        buffer << "</DataArray>" << endl;
    }
}

void VTK::WriteVectorAsVector(stringstream& buffer,
        const Storage3D<dVector3,0>& OutputData, const string VariableName)
{
    // Corresponding PointDataTypes: PDVectors= 2

    buffer << "<DataArray type = \"Float64\" Name = \"" << VariableName <<
                "\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    for(int k = 0; k < OutputData.sizeZ(); ++k)
    for(int j = 0; j < OutputData.sizeY(); ++j)
    for(int i = 0; i < OutputData.sizeX(); ++i)
    {
        buffer << OutputData(i,j,k)[0] << " "
               << OutputData(i,j,k)[1] << " "
               << OutputData(i,j,k)[2] << "\n";
    }
    buffer << "</DataArray>" << endl;
}

void VTK::WriteVectorAsVector(stringstream& buffer,
        const Storage3D<dVector3,1>& OutputData,
        const string VariableName, const int component)
{
    // Corresponding PointDataTypes: PDVectors= 2

        buffer << "<DataArray type = \"Float64\" Name = \""
               << VariableName << "_" << component
               << "\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
        for(int k = 0; k < OutputData.sizeZ(); ++k)
        for(int j = 0; j < OutputData.sizeY(); ++j)
        for(int i = 0; i < OutputData.sizeX(); ++i)
        {
            buffer << OutputData(i,j,k)({component})[0] << " "
                   << OutputData(i,j,k)({component})[1] << " "
                   << OutputData(i,j,k)({component})[2] << "\n";
        }
        buffer << "</DataArray>" << endl;
}

void VTK::WriteSymmetricTensor(stringstream& buffer,
        const Storage3D<dVector6,0>& OutputData, const string VariableName)
{
    // Corresponding PointDataTypes: PDScalars= 1

    vector<long int> dimV {OutputData.sizeX(), OutputData.sizeY(), OutputData.sizeZ()};
    vector<int> compV {0, 1, 2, 3, 4, 5};
    vector<string> compNameV {"xx", "yy", "zz", "yz", "xz", "xy"};

    for (auto it = compV.cbegin(); it != compV.cend(); ++it)
    {
        string compname = "\"" + VariableName + "_" + voigtcon[*it] + "\" ";
        buffer << "<DataArray type = \"Float64\" Name = " << compname <<
                    "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < OutputData.sizeZ(); ++k)
        for(int j = 0; j < OutputData.sizeY(); ++j)
        for(int i = 0; i < OutputData.sizeX(); ++i)
        {
            buffer << OutputData(i,j,k)[*it] << "\n";
        }
        buffer << "</DataArray>" << endl;
    }
}

void VTK::WriteSymmetricTensor(stringstream& buffer,
        const Storage3D<vStrain,0>& OutputData, const string VariableName)
{
    // Corresponding PointDataTypes: PDScalars= 1

    vector<long int> dimV {OutputData.sizeX(), OutputData.sizeY(), OutputData.sizeZ()};
    vector<int> compV {0, 1, 2, 3, 4, 5};
    vector<string> compNameV {"xx", "yy", "zz", "yz", "xz", "xy"};

    for (auto it = compV.cbegin(); it != compV.cend(); ++it)
    {
        string compname = "\"" + VariableName + "_" + voigtcon[*it] + "\" ";
        buffer << "<DataArray type = \"Float64\" Name = " << compname <<
                    "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < OutputData.sizeZ(); ++k)
        for(int j = 0; j < OutputData.sizeY(); ++j)
        for(int i = 0; i < OutputData.sizeX(); ++i)
        {
            buffer << OutputData(i,j,k)[*it] << "\n";
        }
        buffer << "</DataArray>" << endl;
    }
}

void VTK::WriteSymmetricTensor(stringstream& buffer,
        const Storage3D<vStress,0>& OutputData, const string VariableName)
{
    // Corresponding PointDataTypes: PDScalars= 1

    vector<long int> dimV {OutputData.sizeX(), OutputData.sizeY(), OutputData.sizeZ()};
    vector<int> compV {0, 1, 2, 3, 4, 5};
    vector<string> compNameV {"xx", "yy", "zz", "yz", "xz", "xy"};

    for (auto it = compV.cbegin(); it != compV.cend(); ++it)
    {
        string compname = "\"" + VariableName + "_" + voigtcon[*it] + "\" ";
        buffer << "<DataArray type = \"Float64\" Name = " << compname <<
                    "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < OutputData.sizeZ(); ++k)
        for(int j = 0; j < OutputData.sizeY(); ++j)
        for(int i = 0; i < OutputData.sizeX(); ++i)
        {
            buffer << OutputData(i,j,k)[*it] << "\n";
        }
        buffer << "</DataArray>" << endl;
    }
}

void VTK::WriteTensor(stringstream& buffer,
        const Storage3D<int,2>& OutputData, const string VariableName,
        const int componentA, const int componentB)
{
    // Corresponding PointDataTypes: PDScalars= 1
    // Variable name

    buffer << "<DataArray type = \"Float64\" Name = \"" << VariableName <<
                "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
    for(int k = 0; k < OutputData.sizeZ(); ++k)
    for(int j = 0; j < OutputData.sizeY(); ++j)
    for(int i = 0; i < OutputData.sizeX(); ++i)
    {
        buffer << OutputData(i,j,k)({componentA,componentB}) << "\n";
    }
    buffer << "</DataArray>" << endl;
}

void VTK::WriteTensor(stringstream& buffer,
        const Storage3D<double,2>& OutputData, const string VariableName,
        const int componentA, const int componentB)
{
    // Corresponding PointDataTypes: PDScalars= 1

        buffer << "<DataArray type = \"Float64\" Name = \"" << VariableName <<
                    "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < OutputData.sizeZ(); ++k)
        for(int j = 0; j < OutputData.sizeY(); ++j)
        for(int i = 0; i < OutputData.sizeX(); ++i)
        {
            buffer << OutputData(i,j,k)({componentA,componentB}) << "\n";
        }
        buffer << "</DataArray>" << endl;
}

void VTK::WriteQuaternion(stringstream& buffer,
        const Storage3D<Quaternion,0>& OutputData, const string VariableName)
{
    // Corresponding PointDataTypes: PDScalars= 1

    string compname;
    for (int n = 0; n < 4; n++)
    {
        compname = "\"" + VariableName + "_" + to_string(n) + "\" ";
        buffer << "<DataArray type = \"Float64\" Name = " << compname <<
                    "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < OutputData.sizeZ(); ++k)
        for(int j = 0; j < OutputData.sizeY(); ++j)
        for(int i = 0; i < OutputData.sizeX(); ++i)
        {
            buffer << OutputData(i,j,k)[n] << "\n";
        }
        buffer << "</DataArray>" << endl;
    }
}

void VTK::WriteAngles(stringstream& buffer,
        const Storage3D<EulerAngles,0>& OutputData,
        const string VariableName,
        EulerConvention EConvention)
{
    // Corresponding PointDataTypes: PDScalars= 1

    string compname = EulerConventionS[EConvention];
    for (int n = 0; n < 3; n++)
    {
        stringstream converter; 
        converter << compname[n];
        string compname = "\"" + VariableName + "_(" + to_string(n) + ") " + converter.str() + "\" ";
        buffer << "<DataArray type = \"Float64\" Name = " << compname <<
                    "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < OutputData.sizeZ(); ++k)
        for(int j = 0; j < OutputData.sizeY(); ++j)
        for(int i = 0; i < OutputData.sizeX(); ++i)
        {
            buffer << OutputData(i,j,k).Q[n] << "\n";
        }
        buffer << "</DataArray>" << endl;
    }
}

void VTK::WriteMatrix(stringstream& buffer,
        const Storage3D<dMatrix3x3,0>& OutputData, const string VariableName)
{
    // Corresponding PointDataTypes: PDScalars= 1
    string compname;
    for (int n = 0; n < 3; n++)
    for (int m = 0; m < 3; m++)
    {
        compname = "\"" + VariableName + "_" + to_string(n) + to_string(m) + "\" ";
        buffer << "<DataArray type = \"Float64\" Name = " << compname <<
                    "NumberOfComponents=\"1\" format=\"ascii\">" << endl;
        for(int k = 0; k < OutputData.sizeZ(); ++k)
        for(int j = 0; j < OutputData.sizeY(); ++j)
        for(int i = 0; i < OutputData.sizeX(); ++i)
        {
            buffer << OutputData(i,j,k)(n,m) << "\n";
        }
        buffer << "</DataArray>" << endl;
    }
}
}// namespace opensim
