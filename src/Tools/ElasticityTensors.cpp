#include "Tools/Includes.h"
#include "Tools.h"

namespace opensim
{
using namespace std;
vStrain vStrain::Ln() const
{
    dMatrix3x3 M = this->tensor();
    M(0,0) += 1.0;
    M(1,1) += 1.0;
    M(2,2) += 1.0;
    dMatrix3x3 Eigenvectors;
    dMatrix3x3 Eigenvalues;

	Eigenvectors.set_to_zero();
	Eigenvalues.set_to_zero();

	Tools::Eigensystem2(M, Eigenvectors, Eigenvalues);
    //Tools::Eigensystem(M, Eigenvectors, Eigenvalues);
#ifdef DEBUG
    if(fabs(Eigenvalues(0,0)*Eigenvalues(1,1)*Eigenvalues(2,2)) < 1e-6)
    {
        cout << "The Strain's:\n" << this->print() << endl
             << "Matrix:\n"<<M.print()
             << "has at least one zero Eigenvalue:\n"
             << Eigenvalues(0,0)<<endl
             << Eigenvalues(1,1)<<endl
             << Eigenvalues(2,2)<<endl
             << "Calculation of logarithmic strain impossible!\n";
        exit(13);
    }
#endif
    Eigenvalues(0,0) = log(Eigenvalues(0,0));
    Eigenvalues(1,1) = log(Eigenvalues(1,1));
    Eigenvalues(2,2) = log(Eigenvalues(2,2));

    dMatrix3x3 result;
    result.set_to_zero();

    for(int n = 0; n < 3; n++)
    for(int i = 0; i < 3; i++)
    for(int j = i; j < 3; j++)
    {
        result(i,j) += Eigenvalues(n,n)*Eigenvectors(i,n)*Eigenvectors(j,n);
    }

    return result.VoigtStrain();
}

}
