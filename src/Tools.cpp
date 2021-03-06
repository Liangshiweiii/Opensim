#include "Tools.h"
#include "Info.h"
#include "PhaseField.h"
#include "Orientations.h"

namespace opensim
{

using namespace std;
/*
void Tools::WriteAbaqusInput(const PhaseField& Phi, const Orientations& OR)
{
    const int Nx = Phi.Nx;
    const int Ny = Phi.Ny;
    const int Nz = Phi.Nz;
    const double dx = Phi.dx;

    Storage3D<dVector3, 0> Coordinates;                                         // Stores gridpoint coordinates.
    Storage3D<int, 0> ElementNumber;                                            // Stores ElementNumber.
    Storage3D<int, 0> ElementSet;                                               // Stores ElementSet.

    Coordinates.Allocate(Nx+1, Ny+1, Nz+1, 0);
    ElementNumber.Allocate(Nx+2, Ny+2, Nz+2, 0);
    ElementSet.Allocate(Nx+2, Ny+2, Nz+2, 0);

    string GeometryFileName = "geometry.inp";
    string MaterialFileName = "materialinput.inp";
    string GrainFileName = "graindata.inp";

    Info::WriteLine();
    Info::WriteStandard("OpenPhase to Abaqus Parser", "Started");
    Info::WriteStandard("Abaqus geometry file ", GeometryFileName);
    Info::WriteStandard("Abaqus material file ", MaterialFileName);

    // Analyze PhaseField
    int indGrains = 0;
    for(unsigned int n = 0; n < Phi.FieldsStatistics.size(); n++)
    {
        if(Phi.FieldsStatistics[n].Exist)
        {
            indGrains += 1;
        }
    }
    Info::WriteStandard("Individual grains found", std::to_string(indGrains));

    // Determine and write node coordinates
    for (int k =    0; k < Nz+1; k++)
    for (int j =    0; j < Ny+1; j++)
    for (int i =    0; i < Nx+1; i++)
    {
        Coordinates(i,j,k)[0]= i*dx;
        Coordinates(i,j,k)[1]= j*dx;
        Coordinates(i,j,k)[2]= k*dx;
    }

    stringstream outbuffer;

    int elementnumber = 1;
    outbuffer << "*Part, name=PART-1" << endl;
    outbuffer << "*Node, nset=nall" << endl;
    for (int k =    0; k < Nz+1; k++)
    for (int j =    0; j < Ny+1; j++)
    for (int i =    0; i < Nx+1; i++)
    {
        outbuffer << elementnumber << ", "
                  << Coordinates(i,j,k)[0] << ", "
                  << Coordinates(i,j,k)[1] << ", "
                  << Coordinates(i,j,k)[2] << endl;
        elementnumber +=1;
    }

    // Determine and write elements

    Info::WriteStandard("Number of elements (C3D8)", std::to_string(Nx*Ny*Nz));
    outbuffer << "*Element, type=C3D8" << endl;

    elementnumber = 1;
    int elementnodes[8];

    for (int k =    1; k < Nz+1; k++)
    for (int j =    1; j < Ny+1; j++)
    for (int i =    1; i < Nx+1; i++)
    {
        int ii = i;
        int jj = j;
        int kk = k;

        elementnodes[0] = ii + (Nx+1)*(jj-1) + (Nx+1)*(Ny+1)*(kk-1);
        elementnodes[1] = ii + (Nx+1)*(jj-1) + (Nx+1)*(Ny+1)*(kk-1) + 1;
        elementnodes[2] = ii + (Nx+1)*jj     + (Nx+1)*(Ny+1)*(kk-1) + 1;
        elementnodes[3] = ii + (Nx+1)*jj     + (Nx+1)*(Ny+1)*(kk-1);
        elementnodes[4] = ii + (Nx+1)*(jj-1) + (Nx+1)*(Ny+1)*kk;
        elementnodes[5] = ii + (Nx+1)*(jj-1) + (Nx+1)*(Ny+1)*kk + 1;
        elementnodes[6] = ii + (Nx+1)*jj     + (Nx+1)*(Ny+1)*kk + 1;
        elementnodes[7] = ii + (Nx+1)*jj     + (Nx+1)*(Ny+1)*kk;

        ElementNumber(i,j,k) = elementnumber;

        outbuffer << elementnumber << ", ";
        for (int n = 0; n < 7; n++)
        {
            outbuffer << elementnodes[n] << ", ";
        }
        outbuffer << elementnodes[7] << endl;

        // Determine element sets
        // In the interface, the majority phase is used.

        if(Phi.Interface(i-1,j-1,k-1))
        {
            int maxindex = 0.0;
            double maxvalue = 0.0;
            for(auto alpha = Phi.Fields(i-1,j-1,k-1).cbegin();
                    alpha < Phi.Fields(i-1,j-1,k-1).cend(); alpha++)
            {
                if (alpha -> value > maxvalue)
                {
                    maxvalue = alpha -> value;
                    maxindex = alpha -> index;
                }
            }
            ElementSet(i,j,k) = maxindex;
        }
        else
        {
            ElementSet(i,j,k) = Phi.FieldIndex(i-1,j-1,k-1);
        }
        elementnumber += 1;
    }
    outbuffer << "*End Part" << endl;
    outbuffer << "*Assembly, name=Assembly" << endl;
    outbuffer << "*Instance, name=PART-1-1, part=PART-1" << endl;

    // Write element sets
    for (unsigned int ph = 0; ph < Phi.FieldsStatistics.size(); ph++)
    if(Phi.FieldsStatistics[ph].Exist)
    {
        int rowcount = 1;
        outbuffer << "*Elset, elset=Grain" << ph + 1 << endl;
        for (int k =    1; k < Nz+1; k++)
        for (int j =    1; j < Ny+1; j++)
        for (int i =    1; i < Nx+1; i++)
        {
            if((int)ph == ElementSet(i,j,k))
            {
                outbuffer << ElementNumber(i,j,k) << ", ";
                if (rowcount == 15) // Size limitation in Abaqus input files
                {
                    outbuffer << endl;
                    rowcount = 1;
                }
                else
                {
                    rowcount +=1;
                }
            }
        }
        outbuffer << endl << "*Solid Section, elset=Grain" << ph + 1 << ", material=Grain" << ph + 1 << "\n";
    }

    // WriteInterface
    if(Phi.Eta > 0)
    {
        int rowcount = 1;
        outbuffer << "*Elset, elset=InterfaceElements" << endl;
        for (int k =    1; k < Nz+1; k++)
        for (int j =    1; j < Ny+1; j++)
        for (int i =    1; i < Nx+1; i++)
        {
            if(Phi.Interface(i-1,j-1,k-1))
            {
                outbuffer << ElementNumber(i,j,k) << ", ";
                if (rowcount == 15) // Size limitation in Abaqus input files
                {
                    outbuffer << endl;
                    rowcount = 1;
                }
                else
                {
                    rowcount +=1;
                }
            }
        }
        outbuffer << endl;
        for (unsigned int ph = 0; ph < Phi.FieldsStatistics.size(); ph++)
        if(Phi.FieldsStatistics[ph].Exist)
        {
            int rowcount = 1;
            outbuffer << "*Elset, elset=Interface_" << ph << endl;
            for (int k =    1; k < Nz+1; k++)
            for (int j =    1; j < Ny+1; j++)
            for (int i =    1; i < Nx+1; i++)
            {
                if(Phi.Interface(i-1,j-1,k-1) and Phi.Fields(i-1,j-1,k-1).get(ph) >= 0.5)
                {
                    outbuffer << ElementNumber(i,j,k) << ", ";
                    if (rowcount == 15) // Size limitation in Abaqus input files
                    {
                        outbuffer << endl;
                        rowcount = 1;
                    }
                    else
                    {
                        rowcount +=1;
                    }
                }
            }
            outbuffer << endl;
        }
    } // end interface existing

    int rowcount = 1;
    outbuffer << "*Elset, elset=MIDDLEELEMENTS" << endl;
    for (int i =    1; i < Nz+1; i++)
    {
        outbuffer << ElementNumber(Nx/2,Ny/2+1, i) << ", ";
        if (rowcount == 15) // Size limitation in Abaqus input files
        {
            outbuffer << endl;
            rowcount = 1;
        }
        else
        {
            rowcount +=1;
        }
    }
    outbuffer << endl;

    outbuffer << "*End Instance" << endl;
    outbuffer << "*End Assembly" << endl;
    outbuffer << "*Include, input=materialinput.inp" << endl;
    //outbuffer << "*Include, input=boundary.inp" << endl;

    ofstream Geometry_file(GeometryFileName.c_str());
    Geometry_file << outbuffer.rdbuf();
    Geometry_file.close();

    outbuffer.str(std::string());                                               // wipe
    for (unsigned int ph = 0; ph < Phi.FieldsStatistics.size(); ph++)
    if(Phi.FieldsStatistics[ph].Exist)
    {
        outbuffer << "*Material, name=Grain" << ph + 1 << endl
                  << "*Depvar" << endl
                  << "     53" << endl
                  << "*User Material, constants = 2" << endl
                  << ph + 1  << " , 2" << endl;
    }

    ofstream Material_file(MaterialFileName.c_str());
    Material_file << outbuffer.rdbuf();
    Material_file.close();

    Info::WriteStandard("Abaqus grain data file ", GrainFileName);

    // Writes seperate file with Euler angles
    ofstream Grain_file(GrainFileName.c_str());
    outbuffer.str(std::string());                                               // wipe

    for (unsigned int ph = 0; ph < Phi.FieldsStatistics.size(); ph++)
    if(Phi.FieldsStatistics[ph].Exist)
    {
        Angles tempOrient = OR.GrainEulerAngles[ph].get_degree();
        outbuffer << "Grain : " << ph + 1 << " : " <<
                           tempOrient.Q[0] <<  " : " <<
                           tempOrient.Q[1] <<  " : " <<
                           tempOrient.Q[2] << endl;
    }

    Grain_file << outbuffer.rdbuf();
    Grain_file.close();

    Info::WriteLine();
}
*/
void Tools::Eigensystem(dMatrix3x3& M, dMatrix3x3& Eigenvectors, dMatrix3x3& Eigenvalues)
{
#ifdef DEBUG
    ///M MUST BE _SYMMETRIC_
    bool symmetric = true;
    for(int i = 0; i < 3; i++)
    for(int j = i; j < 3; j++)
    {
        if(M(i,j) != M(j,i))
            symmetric = false;
    }
    if(not symmetric)
    {
        Info::WriteExit("M:\n"
                        + M.print() + "is not symmetric!", "Tools", "Eigensystem()");
        exit(13);
    }
#endif
    int max_iterations=1000000;
    double Check = 1e3, delta=1e-4;
    dMatrix3x3 tmpMatrix, M2,M4,M8,
    //M16,M32,M64,M128, M256, M512,
    MFinal;
    M2=M*M;
    M4=M2*M2;
    M8=M4*M4;
    //M16=M8*M8;
    //M32=M16*M16;
    //M64=M32*M32;
    //M128=M64*M64;
    //M256=M128*M128;
    //M512=M256*M256;
    MFinal=M;

    dVector3 u0;
    u0.set_to_unitX();
    dVector3 tmp;
    tmp = M*u0;                                                                 //tmp = M.u0
    tmp.normalize();
    dVector3 u = tmp;

    dVector3 v = u-u0;                                                          //{u[0] - u0[0], u[1] - u0[1], u[2] - u0[2]};
    if(v.length() < 1e-3)
        v.setY(1.0);
    v.normalize();

    int n = 1;
    while(Check >= delta && fabs(2.0 - Check) > delta)
    {
        tmp.set_to_zero();
        tmp -= u;

        u = MFinal*u;                                                           //u_i=M^32 . u_(i-1), high power of M accelerates convergency
        u.normalize();

        tmp += u;

        Check = tmp.length();                                                   //TODO: to avoid sqrt() calculation it's possible to use smth. like a Lengthsqr() fct!

        n++;
        if(n >= max_iterations and MFinal != M8)
        {
            n=0;
            MFinal = M8;
        }

        if(n >= max_iterations and MFinal == M8)
        {
            stringstream msg;
            msg << "Eigensystem() could not converge! M:\n " << M.print()
                << "Last suggested solution for the first eigenvector: "
                << u.print() << "; last correction: " << Check <<" > "
                << delta << endl;
            Info::WriteExit( msg.str(), "Tools", "Eigensystem()");
            exit(13);
        }
    }
//    cout<<n-1<<" iterations, Check: "<<Check<<endl;

    double v_cdot_u = u*v;
    v -= (u*v_cdot_u);
    v.normalize();

    n = 1;
    Check = 1e3;
    MFinal = M;
    while(Check >= delta && fabs(2.0 - Check) > delta)
    {
        tmp.set_to_zero();
        tmp -= v;
        v = MFinal*v;                                                           //v_i=M^32 . v_(i-1) ; it's better not to use to high powers of M here (-;

        v_cdot_u = u*v;
        v -= (u*v_cdot_u);
        v.normalize();

        tmp += v;

        Check = tmp.length();                                                   //TODO: to avoid sqrt() calculation it's possible to use smth. like a Lengthsqr() fct!

        n++;
        if(n >= max_iterations and MFinal != M8)
        {
            n=0;
            MFinal = M8;
        }
        if(n >= max_iterations and MFinal == M8)
        {
            stringstream msg;
            msg << "Eigensystem() could not converge! M:\n " << M.print()
                << "Last suggested solution for the second eigenvector: "
                << v.print() << "; last correction: " << Check <<" > "
                << delta << endl;
            Info::WriteExit(msg.str(), "Tools", "Eigensystem()");
            exit(13);
        }
    }
//    cout<<n-1<<" iterations, Check: "<<Check<<endl;

    dVector3 w;
    w = u.cross(v);

    dMatrix3x3 Eigenvectors_tr;                                                 //Row Eigenvectors
    Eigenvectors_tr(0,0) = u.getX();
    Eigenvectors_tr(0,1) = u.getY();
    Eigenvectors_tr(0,2) = u.getZ();
    Eigenvectors_tr(1,0) = v.getX();
    Eigenvectors_tr(1,1) = v.getY();
    Eigenvectors_tr(1,2) = v.getZ();
    Eigenvectors_tr(2,0) = w.getX();
    Eigenvectors_tr(2,1) = w.getY();
    Eigenvectors_tr(2,2) = w.getZ();

    Eigenvectors = Eigenvectors_tr.transposed();                                //Column Eigenvectors

    tmpMatrix = Eigenvectors_tr * M;
    Eigenvalues = tmpMatrix*Eigenvectors;
    for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
        if(i != j)
            Eigenvalues(i,j) = 0.0;
    }

//    if(Eigenvectors(0,0)<1.0)
//    {
//        cout<<n-1<<" iterations, Check: "<<Check<<endl;
//        cout<<"Eigenvalues:\n"<<Eigenvalues.print()<<endl;
//        cout<<"Eigenvectors:\n"<<Eigenvectors.print()<<endl;
//    }
}

void Tools::Eigensystem2(dMatrix3x3& M, dMatrix3x3& Eigenvectors, dMatrix3x3& Eigenvalues, bool calculateEigenVectors)
{
	if (M(0, 1) != M(1, 0) or M(0, 2) != M(2, 0) or M(1, 2) != M(2, 1))
    {
        stringstream msg;
        msg << "Matrix is not symmetric M:\n " << M.print() << endl;
        Info::WriteWarning(msg.str(), "Tools", "Eigensystem2()");
    }
	double p1 = pow(M(0, 1), 2) + pow(M(0, 2), 2) + pow(M(1, 2), 2);
	if (p1 == 0)
	{
		Eigenvalues(0, 0) = M(0, 0);
		Eigenvalues(1, 1) = M(1, 1);
		Eigenvalues(2, 2) = M(2, 2);
	}
	else
	{
		double q = M.trace() / 3.0;
		double p2 = pow((M(0, 0) - q), 2) + pow((M(1, 1) - q), 2) + pow((M(2, 2) - q), 2) + 2 * p1;
		double p = sqrt(p2 / 6.0);

		if (p != 0.0)
		{
			dMatrix3x3 I;
			I.set_to_unity();
			dMatrix3x3 B = (M - I * q) * (1 / p);
			double r = 0.5 * B.determinant();

			double phi = 0;
			if (r <= -1.0)
				phi = Pi / 3.0;
			else if (r < 1.0)
				phi = acos(r) / 3.0;

			//the eigenvalues satisfy eig3 <= eig2 <= eig1
			Eigenvalues(0, 0) = q + 2.0 * p * cos(phi);
			Eigenvalues(2, 2) = q + 2.0 * p * cos(phi + (2.0*Pi / 3.0));
			Eigenvalues(1, 1) = 3.0 * q - Eigenvalues(0, 0) - Eigenvalues(2, 2);
		}
		else
		{
			Eigenvalues(0, 0) = 1;
			Eigenvalues(2, 2) = 1;
			Eigenvalues(1, 1) = 1;
		}
	}

	Eigenvectors.set_to_unity();
	if (calculateEigenVectors)
	{
		if (Eigenvalues(0, 0) != Eigenvalues(1, 1) and Eigenvalues(0, 0) != Eigenvalues(2, 2))
			for (int n = 0; n < 3; n++)
			{
				dVector3 row0({ M(0,0) - Eigenvalues(n,n), M(0,1), M(0,2) });
				dVector3 row1({ M(0,1), M(1,1) - Eigenvalues(n,n), M(0,2) });
				dVector3 row2({ M(0,2), M(1,2), M(2,2) - Eigenvalues(n,n) });

				dVector3 r0xr1 = row0.cross(row1);
				dVector3 r0xr2 = row0.cross(row2);
				dVector3 r1xr2 = row1.cross(row2);

				double d0 = r0xr1 * r0xr1;
				double d1 = r0xr2 * r0xr2;
				double d2 = r1xr2 * r1xr2;
				double dMax = d0;
				int iMax = 0;

				if (d1 > dMax)
				{
					dMax = d1;
					iMax = 1;
				}
				if (d2 > dMax)
				{
					iMax = 2;
				}
				if (iMax == 0 and d0 > 0.0)
				{
					dVector3 temp = r0xr1 / sqrt(d0);
					Eigenvectors(n, 0) = temp[0];
					Eigenvectors(n, 1) = temp[1];
					Eigenvectors(n, 2) = temp[2];
				}

				else if (iMax == 1 and d1 > 0.0)
				{
					dVector3 temp = r0xr2 / sqrt(d1);
					Eigenvectors(n, 0) = temp[0];
					Eigenvectors(n, 1) = temp[1];
					Eigenvectors(n, 2) = temp[2];
				}
				else if (d2 > 0.0)
				{
					dVector3 temp = r1xr2 / sqrt(d2);
					Eigenvectors(n, 0) = temp[0];
					Eigenvectors(n, 1) = temp[1];
					Eigenvectors(n, 2) = temp[2];
				}
			}
	}

			//cout << "Eigenvalues\n" << Eigenvalues.print() << "\nEigenvectors\n" << Eigenvectors.print() << "\n" << endl;
}

void Tools::Decompose(dMatrix3x3& Deformation, dMatrix3x3& Rot, dMatrix3x3& Strain)
{
    dMatrix3x3 Deformation_tr, M, Eigenvectors, Eigenvectors_tr,
               Eigenvalues,tmp, Strain_inverse;
    Deformation_tr = Deformation.transposed();
    M = Deformation_tr * Deformation;
    Eigensystem(M,Eigenvectors,Eigenvalues);

    for(int i=0; i<3; i++)
    {
        Eigenvalues(i,i) = sqrt(Eigenvalues(i,i));
    }
    Eigenvectors_tr = Eigenvectors.transposed();
    tmp = Eigenvalues * Eigenvectors_tr;
    Strain = Eigenvectors*tmp;                                                  //Strain = Eigenvectors_column . sqrt(Eigenvalues) . Eigenvectors_row
    Strain_inverse = Strain.inverted();
    Rot = Deformation*Strain_inverse;                                           //Rotate= InverseStrain . Deformation
}

}// namespace opensim
