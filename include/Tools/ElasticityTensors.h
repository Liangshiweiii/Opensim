#ifndef ELASTICITYTENSORS_H
#define ELASTICITYTENSORS_H

//#include "Base/Includes.h"
#include <cfloat>
#include <cmath>
#include <cstring>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "Tools/ElasticityTensors.h"

namespace opensim
{

class dVector3;
class dMatrix3x3;
class dVector6;
class dMatrix6x6;
class vStrain;
class vStress;

class dMatrix3x3
{
 public:

    dMatrix3x3(){};

    dMatrix3x3(const dMatrix3x3& rhs)
    {
        memmove(storage, rhs.const_data(), 9*sizeof(double));
    };

    dMatrix3x3(std::initializer_list<double> matinit)
    {
        storage[0][0] = 0.0;
        storage[1][0] = 0.0;
        storage[2][0] = 0.0;
        storage[0][1] = 0.0;
        storage[1][1] = 0.0;
        storage[2][1] = 0.0;
        storage[0][2] = 0.0;
        storage[1][2] = 0.0;
        storage[2][2] = 0.0;

#ifdef DEBUG
        if (matinit.size() != 9)
        {
            std::cout << "Error in dMatrix3x3::constructor()\n"
                      << "Initialization list size beyond storage range."
                      << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        int ii = 0;
        int jj = 0;
        for (auto it = matinit.begin(); it != matinit.end(); it++)
        {
            if (ii == 3)
            {
                ii =  0;
                jj += 1;
            }
            storage[jj][ii] = *it;
            ii += 1;
        }
    }

    double& operator()(const int i, const int j)
    {
#ifdef DEBUG
        if(i > 2 or j > 2)
        {
            std::cout << "Error in dMatrix3x3::operator()\n"
                 << "Access beyond storage range. (i, j) = "
                 << i <<", "<< j << " > (2, 2)"
                 << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i][j];
    };
    double const& operator()(const int i, const int j) const
    {
#ifdef DEBUG
        if(i > 2 or j > 2)
        {
            std::cout << "Error in dMatrix3x3::operator()\n"
                 << "Access beyond storage range. (i, j) = "
                 << i <<", "<< j << " > (2, 2)"
                 << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i][j];
    };
    void set_to_zero(void)
    {
        memset(storage, 0, 9*sizeof(double));
    };
    void set_to_unity(void)
    {
        memset(storage, 0, 9*sizeof(double));
        storage[0][0] = 1.0;
        storage[1][1] = 1.0;
        storage[2][2] = 1.0;
    };
    dMatrix3x3& operator=(const dMatrix3x3& rhs)
    {
        memmove(storage, rhs.const_data(), 9*sizeof(double));
        return *this;
    };

    bool operator==(const dMatrix3x3& rhs)
    {
        for(int i=0; i < 3; i++)
        for(int j=0; j < 3; j++)
        {
            if(storage[i][j]!=rhs(i,j))
                return false;
        }
        return true;
    };

    bool operator!=(const dMatrix3x3& rhs)
    {
        for(int i=0; i < 3; i++)
        for(int j=0; j < 3; j++)
        {
            if(storage[i][j]!=rhs(i,j))
                return true;
        }
        return false;
    };

    dMatrix3x3 operator*(const double m) const
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) = storage[i][j]*m;
        }
        return tmp;
    };
    dMatrix3x3 operator/(const double m) const
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) = storage[i][j]/m;
        }
        return tmp;
    };
    dVector3 operator*(const dVector3& rhs) const;
    // the method is implemented ouside of the class (see below)
   /* {
        dVector3 tmp;
        memset(tmp, 0, 3*sizeof(double));
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp[i] += storage[i][j]*rhs[j];
        }
        return tmp;
    }; */
    dMatrix3x3 operator*(const dMatrix3x3& rhs) const
    {
        dMatrix3x3 tmp;
        tmp.set_to_zero();
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        for(int k = 0; k < 3; k++)
        {
            tmp(i,j) += storage[i][k]*rhs(k,j);
        }
        return tmp;
    };
    dMatrix3x3 operator+(const dMatrix3x3& rhs) const
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) = storage[i][j] + rhs(i,j);
        }
        return tmp;
    };
    dMatrix3x3 operator-(const dMatrix3x3& rhs) const
    {
        dMatrix3x3 tmp;
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp(i,j) = storage[i][j] - rhs(i,j);
        }
        return tmp;
    };

    dMatrix6x6 outer(const dMatrix3x3& rhs) const;

    dMatrix3x3& operator+=(const dMatrix3x3& rhs)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            storage[i][j] += rhs(i,j);
        }
        return *this;
    };
    dMatrix3x3& operator-=(const dMatrix3x3& rhs)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            storage[i][j] -= rhs(i,j);
        }
        return *this;
    };
    dMatrix3x3& operator*=(const double m)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            storage[i][j] *= m;
        }
        return *this;
    };
    dMatrix3x3& operator/=(const double m)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            storage[i][j] /= m;
        }
        return *this;
    };
    inline double determinant(void) const
    {
        return (storage[0][0]*storage[1][1]*storage[2][2] +
                storage[0][1]*storage[1][2]*storage[2][0] +
                storage[0][2]*storage[1][0]*storage[2][1] -
                storage[0][2]*storage[1][1]*storage[2][0] -
                storage[0][1]*storage[1][0]*storage[2][2] -
                storage[0][0]*storage[1][2]*storage[2][1]);
    };
    dMatrix3x3& invert(void)
    {
        dMatrix3x3 tmp;

        double detInv = determinant();

        if(detInv != 0.0) detInv = 1.0/detInv;
        else
        {
            std::cout << "dMatrix3x3:\n" << this->print() << "Is not Invertible!"
                      << std::endl;
            std::cout << tmp.print() << std::endl;
            exit(13);
        }

        tmp(0,0) = (storage[1][1]*storage[2][2] -
                    storage[1][2]*storage[2][1])*detInv;
        tmp(1,0) =-(storage[1][0]*storage[2][2] -
                    storage[1][2]*storage[2][0])*detInv;
        tmp(2,0) = (storage[1][0]*storage[2][1] -
                    storage[1][1]*storage[2][0])*detInv;
        tmp(0,1) =-(storage[0][1]*storage[2][2] -
                    storage[0][2]*storage[2][1])*detInv;
        tmp(1,1) = (storage[0][0]*storage[2][2] -
                    storage[0][2]*storage[2][0])*detInv;
        tmp(2,1) =-(storage[0][0]*storage[2][1] -
                    storage[0][1]*storage[2][0])*detInv;
        tmp(0,2) = (storage[0][1]*storage[1][2] -
                    storage[1][1]*storage[0][2])*detInv;
        tmp(1,2) =-(storage[0][0]*storage[1][2] -
                    storage[0][2]*storage[1][0])*detInv;
        tmp(2,2) = (storage[0][0]*storage[1][1] -
                    storage[0][1]*storage[1][0])*detInv;

        memmove(storage, tmp.data(), 9*sizeof(double));
        return *this;
    };
    dMatrix3x3 inverted(void) const
    {
        dMatrix3x3 tmp;

        double detInv = determinant();

        if(detInv != 0.0) detInv = 1.0/detInv;
        else
        {
            std::cout << "dMatrix3x3:\n" << this->print() << " Is not Invertible!"
                      << std::endl;
            std::cout << tmp.print() << std::endl;
            exit(13);
        }

        tmp(0,0) = (storage[1][1]*storage[2][2] -
                    storage[1][2]*storage[2][1])*detInv;
        tmp(1,0) =-(storage[1][0]*storage[2][2] -
                    storage[1][2]*storage[2][0])*detInv;
        tmp(2,0) = (storage[1][0]*storage[2][1] -
                    storage[1][1]*storage[2][0])*detInv;
        tmp(0,1) =-(storage[0][1]*storage[2][2] -
                    storage[0][2]*storage[2][1])*detInv;
        tmp(1,1) = (storage[0][0]*storage[2][2] -
                    storage[0][2]*storage[2][0])*detInv;
        tmp(2,1) =-(storage[0][0]*storage[2][1] -
                    storage[0][1]*storage[2][0])*detInv;
        tmp(0,2) = (storage[0][1]*storage[1][2] -
                    storage[1][1]*storage[0][2])*detInv;
        tmp(1,2) =-(storage[0][0]*storage[1][2] -
                    storage[0][2]*storage[1][0])*detInv;
        tmp(2,2) = (storage[0][0]*storage[1][1] -
                    storage[0][1]*storage[1][0])*detInv;

        return tmp;
    };
    dMatrix3x3& transpose(void)
    {
        double tmp[3][3];
        memset(tmp, 0, 9*sizeof(double));

        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            tmp[i][j] = storage[j][i];
        }
        memmove(storage, tmp, 9*sizeof(double));
        return *this;
    };
    dMatrix3x3 transposed(void) const
    {
        dMatrix3x3 TempMat;

        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            TempMat(i,j) = storage[j][i];
        }
        return TempMat;
    };
    dMatrix3x3& rotate(const dMatrix3x3& RotationMatrix)
    {
        dMatrix3x3 TempMat;
        //void const* const rightptr = storage;
        //void * leftptr = TempMat.data();

        memmove(TempMat.data(), storage, 9*sizeof(double));

        TempMat = RotationMatrix * (TempMat * RotationMatrix.transposed());

        memmove(storage, TempMat.data(), 9*sizeof(double));
        return *this;
    };
    dMatrix3x3 rotated(const dMatrix3x3& RotationMatrix) const
    {
        dMatrix3x3 Out;
        memmove(Out.data(), storage, 9*sizeof(double));

        Out = RotationMatrix * (Out * RotationMatrix.transposed());

        return Out;
    };
    double double_contract(const dMatrix3x3& rHS) const
    {
        double tmp = 0.0;
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            tmp += storage[i][j]*rHS(i,j);
        }
        return tmp;
    };
    double trace(void) const
    {
        return storage[0][0] + storage[1][1] + storage[2][2];
    };
    dMatrix3x3 getsym(void) const
    {
        dMatrix3x3 TempMat;

        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            TempMat(i,j) = storage[i][j] + storage[j][i];
        }
        TempMat *= 0.5;

        return TempMat;
    };
    dMatrix3x3 getskew(void) const
    {
        dMatrix3x3 TempMat;

        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            TempMat(i,j) = storage[i][j] - storage[j][i];
        }
        TempMat *= 0.5;

        return TempMat;
    };
    double norm(void) const  /// Frobenius norm
    {
        double tmp = 0.0;
        for(int i = 0; i < 3; i++)
        {
            tmp += storage[i][0] * storage[i][0]
                 + storage[i][1] * storage[i][1]
                 + storage[i][2] * storage[i][2];
        }

        return sqrt(tmp);
    };
    std::string print(void) const
    {
        std::stringstream out;
        for(int i = 0; i < 3; i++)
        {
            out << "||" << std::setprecision(6) << std::right
                        << std::setw(8) << storage[i][0] << " "
                        << std::setw(8) << storage[i][1] << " "
                        << std::setw(8) << storage[i][2] << "||\n";
        }
        return out.str();
    };

    void read(std::fstream& inp)
    {
        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        {
            inp >> storage[i][j];
        }
    };
    double * data(void)
    {
        return &storage[0][0];
    };
    double const * const_data(void) const
    {
        return &storage[0][0];
    };
    inline dVector6 VoigtVector() const;
    inline vStrain  VoigtStrain() const;
    inline vStress  VoigtStress() const;
 protected:
 private:

    double storage[3][3];
};

class dVector3
{
 public:

    dVector3(){};

    dVector3(std::initializer_list<double> vecinit)
    {
        storage[0] = 0.0;
        storage[1] = 0.0;
        storage[2] = 0.0;

#ifdef DEBUG
        if (vecinit.size() != 3)
        {
            std::cout << "Error in dVector3::constructor()\n"
                      << "Initialization list size beyond storage range."
                      << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        int ii = 0;
        for (auto it = vecinit.begin(); it != vecinit.end(); it++)
        {
            storage[ii] = *it;
            ii += 1;
        }
    }

    double& operator[](const int i)
    {
#ifdef DEBUG
        if(i > 2)
        {
            std::cout << "Error in dVector3::operator[]\n"
                      << "Access beyond storage range. i = "
                      << i << " > 2"
                      << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i];
    };
    double const& operator[](const int i) const
    {
#ifdef DEBUG
        if(i > 2)
        {
            std::cout << "Error in dVector3::operator[]\n"
                      << "Access beyond storage range. i = "
                      << i << " > 2"
                      << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i];
    };

    double getX(void) const
    {
        return storage[0];
    };

    void setX(const double newX)
    {
        storage[0] = newX;
    };

    double getY(void) const
    {
        return storage[1];
    };

    void setY(const double newY)
    {
        storage[1] = newY;
    };

    double getZ(void) const
    {
        return storage[2];
    };

    void setZ(const double newX)
    {
        storage[2] = newX;
    };

    void set_to_zero(void)
    {
        memset(storage, 0, 3*sizeof(double));
    };
    void set_to_unitX(void)
    {
        storage[0] = 1.0;
        storage[1] = 0.0;
        storage[2] = 0.0;
    };
    void set_to_unitY(void)
    {
        storage[0] = 0.0;
        storage[1] = 1.0;
        storage[2] = 0.0;
    };
    void set_to_unitZ(void)
    {
        storage[0] = 0.0;
        storage[1] = 0.0;
        storage[2] = 1.0;
    };
    dVector3 operator*(const double m) const
    {
        dVector3 tmp;
        tmp[0] = storage[0]*m;
        tmp[1] = storage[1]*m;
        tmp[2] = storage[2]*m;
        return tmp;
    };
    dVector3 operator/(const double m) const
    {
        dVector3 tmp;
        tmp[0] = storage[0]/m;
        tmp[1] = storage[1]/m;
        tmp[2] = storage[2]/m;
        return tmp;
    };
    double operator*(const dVector3& rhs) const
    {
        return storage[0]*rhs[0] + storage[1]*rhs[1] + storage[2]*rhs[2];
    };
    double abs() const
    {
        return sqrt(storage[0]*storage[0] +
                    storage[1]*storage[1] +
                    storage[2]*storage[2]);
    };
    dVector3 cross(const dVector3& rhs) const
    {
        dVector3 tmp;
        tmp[0] = storage[1] * rhs[2] - storage[2] * rhs[1];
        tmp[1] = storage[2] * rhs[0] - storage[0] * rhs[2];
        tmp[2] = storage[0] * rhs[1] - storage[1] * rhs[0];
        return tmp;
    };
    dVector3 operator+(const dVector3& rhs) const
    {
        dVector3 tmp;
        tmp[0] = storage[0] + rhs[0];
        tmp[1] = storage[1] + rhs[1];
        tmp[2] = storage[2] + rhs[2];
        return tmp;
    };
    dVector3 operator-(const dVector3& rhs) const
    {
        dVector3 tmp;
        tmp[0] = storage[0] - rhs[0];
        tmp[1] = storage[1] - rhs[1];
        tmp[2] = storage[2] - rhs[2];
        return tmp;
    };
    dVector3& operator*=(const double m)
    {
        storage[0] *= m;
        storage[1] *= m;
        storage[2] *= m;
        return *this;
    };
    dVector3& operator/=(const double m)
    {
        storage[0] /= m;
        storage[1] /= m;
        storage[2] /= m;
        return *this;
    };
    dVector3& operator-=(const dVector3& rhs)
    {
        storage[0] = storage[0] - rhs[0];
        storage[1] = storage[1] - rhs[1];
        storage[2] = storage[2] - rhs[2];
        return *this;
    };
    dVector3& operator+=(const dVector3& rhs)
    {
        storage[0] = storage[0] + rhs[0];
        storage[1] = storage[1] + rhs[1];
        storage[2] = storage[2] + rhs[2];
        return *this;
    };
    dVector3& operator=(const dVector3& rhs)
    {
        memmove(data(), rhs.const_data(), 3*sizeof(double));
        return *this;
    };
    dVector3& operator=(const double rhs[3])
    {
        memmove(storage, rhs, 3*sizeof(double));
        return *this;
    };
    double length(void) const
    {
        return sqrt(storage[0]*storage[0] +
                    storage[1]*storage[1] +
                    storage[2]*storage[2]);
    };
    dVector3& normalize(void)
    {
        double norm = sqrt(storage[0]*storage[0] +
                           storage[1]*storage[1] +
                           storage[2]*storage[2]);
        if(norm != 0.0)
        {
            double norm_inv = 1.0/norm;
            storage[0] *= norm_inv;
            storage[1] *= norm_inv;
            storage[2] *= norm_inv;
        }
        return *this;
    };
    dVector3 normalized(void) const
    {
        double norm = sqrt(storage[0]*storage[0] +
                           storage[1]*storage[1] +
                           storage[2]*storage[2]);
        dVector3 tmp;
        if(norm != 0.0)
        {
            double norm_inv = 1.0/norm;
            tmp[0] = storage[0]*norm_inv;
            tmp[1] = storage[1]*norm_inv;
            tmp[2] = storage[2]*norm_inv;
        }
        else
        {
            tmp.set_to_zero();
        }
        return tmp;
    };
    dVector3& rotate(const dMatrix3x3& RotationMatrix)
    {
        dVector3 Out;
        Out.set_to_zero();
        for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
        {
            Out[i] += RotationMatrix(i,j) * storage[j];
        }
        memmove(storage, Out.data(), 3*sizeof(double));
        return *this;
    };
    dVector3 rotated(const dMatrix3x3& RotationMatrix) const
    {
        dVector3 Out;
        Out.set_to_zero();
        for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
        {
            Out[i] += RotationMatrix(i,j) * storage[j];
        }
        return Out;
    };
    dVector3 Xreflected(void) const
    {
        dVector3 Out;
        for(int i = 0; i < 3; ++i)
        {
            Out[i] = storage[i];
        }
        Out[0] *= -1.0;
        return Out;
    };
    dVector3 Yreflected(void) const
    {
        dVector3 Out;
        for(int i = 0; i < 3; ++i)
        {
            Out[i] = storage[i];
        }
        Out[1] *= -1.0;
        return Out;
    };
    dVector3 Zreflected(void) const
    {
        dVector3 Out;
        for(int i = 0; i < 3; ++i)
        {
            Out[i] = storage[i];
        }
        Out[2] *= -1.0;
        return Out;
    };
    std::string print(void) const
    {
        std::stringstream out;

        out << "(" << storage[0] << ", "
                   << storage[1] << ", "
                   << storage[2] << ")";
        return out.str();
    };
    double* data(void)
    {
        return storage;
    };
    const double* const_data(void) const
    {
        return storage;
    };
 protected:
 private:

    double storage[3];
};

inline dVector3 dMatrix3x3::operator*(const dVector3& rhs) const
{
    dVector3 tmp;
    tmp.set_to_zero();
    for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
    {
        tmp[i] += storage[i][j]*rhs[j];
    }
    return tmp;
}

class dVector6
{
 public:

    dVector6(){};

    dVector6(std::initializer_list<double> vecinit)
    {
#ifdef DEBUG
        if (vecinit.size() != 6)
        {
            std::cout << "Error in dVector6::constructor()\n"
                      << "Initialization list size beyond storage range."
                      << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        int ii = 0;
        for (auto it = vecinit.begin(); it != vecinit.end(); it++)
        {
            storage[ii] = *it;
            ii += 1;
        }
    }

    double& operator[](const int i)
    {
#ifdef DEBUG
        if(i > 5)
        {
            std::cout << "Error in dVector6::operator[]\n"
                      << "Access beyond storage range. i = "
                      << i << " > 5"
                      << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i];
    };
    double const& operator[](const int i) const
    {
#ifdef DEBUG
        if(i > 5)
        {
            std::cout << "Error in dVector6::operator[]\n"
                      << "Access beyond storage range. i = "
                      << i << " > 5"
                      << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i];
    };
    void set_to_zero(void)
    {
        memset(storage, 0, 6*sizeof(double));
    };
    void set_to_unity(void)
    {
        storage[0] = 1.0;
        storage[1] = 1.0;
        storage[2] = 1.0;
        storage[3] = 0.0;
        storage[4] = 0.0;
        storage[5] = 0.0;
    };
    double trace(void) const
    {
        return storage[0] + storage[1] + storage[2];
    };
    double max_abs(void) const
    {
        // Returns maximum absolute value
        double tempmax = 0.0;
        for (int i = 0; i < 6; i++)
        {
            if (abs(storage[i]) > tempmax)
            {
                tempmax = abs(storage[i]);
            }
        }
        return tempmax;
    };
//    double norm(void) const  /// Frobenius norm
//    {
//        double tmp = storage[0] * storage[0]
//                   + storage[1] * storage[1]
//                   + storage[2] * storage[2]
//                   + storage[3] * storage[3] * 2.0
//                   + storage[4] * storage[4] * 2.0
//                   + storage[5] * storage[5] * 2.0;
//
//        return sqrt(tmp);
//    };
    dVector6 operator*(const double m) const
    {
        dVector6 tmp;
        tmp[0] = storage[0]*m;
        tmp[1] = storage[1]*m;
        tmp[2] = storage[2]*m;
        tmp[3] = storage[3]*m;
        tmp[4] = storage[4]*m;
        tmp[5] = storage[5]*m;
        return tmp;
    };
    dVector6 operator/(const double m) const
    {
        dVector6 tmp;
        tmp[0] = storage[0]/m;
        tmp[1] = storage[1]/m;
        tmp[2] = storage[2]/m;
        tmp[3] = storage[3]/m;
        tmp[4] = storage[4]/m;
        tmp[5] = storage[5]/m;
        return tmp;
    };
    dVector6& operator*=(const double m)
    {
        storage[0] *= m;
        storage[1] *= m;
        storage[2] *= m;
        storage[3] *= m;
        storage[4] *= m;
        storage[5] *= m;
        return *this;
    };

    double operator*(const dVector6& rhs)
    {
        double tmp = 0.0;
        for(int i = 0; i < 6; i++)
        {
            tmp += storage[i]*rhs[i];
        }
        return tmp;
    };

    dVector6& operator/=(const double m)
    {
        storage[0] /= m;
        storage[1] /= m;
        storage[2] /= m;
        storage[3] /= m;
        storage[4] /= m;
        storage[5] /= m;

        return *this;
    };
    dVector6 operator+(const dVector6& rhs) const
    {
        dVector6 tmp;
        tmp[0] = storage[0] + rhs[0];
        tmp[1] = storage[1] + rhs[1];
        tmp[2] = storage[2] + rhs[2];
        tmp[3] = storage[3] + rhs[3];
        tmp[4] = storage[4] + rhs[4];
        tmp[5] = storage[5] + rhs[5];
        return tmp;
    };
    dVector6& operator+=(const dVector6& rhs)
    {
        storage[0] = storage[0] + rhs[0];
        storage[1] = storage[1] + rhs[1];
        storage[2] = storage[2] + rhs[2];
        storage[3] = storage[3] + rhs[3];
        storage[4] = storage[4] + rhs[4];
        storage[5] = storage[5] + rhs[5];
        return *this;
    };
    dVector6 operator-(const dVector6& rhs) const
    {
        dVector6 tmp;
        tmp[0] = storage[0] - rhs[0];
        tmp[1] = storage[1] - rhs[1];
        tmp[2] = storage[2] - rhs[2];
        tmp[3] = storage[3] - rhs[3];
        tmp[4] = storage[4] - rhs[4];
        tmp[5] = storage[5] - rhs[5];
        return tmp;
    };
    dVector6& operator-=(const dVector6& rhs)
    {
        storage[0] = storage[0] - rhs[0];
        storage[1] = storage[1] - rhs[1];
        storage[2] = storage[2] - rhs[2];
        storage[3] = storage[3] - rhs[3];
        storage[4] = storage[4] - rhs[4];
        storage[5] = storage[5] - rhs[5];
        return *this;
    };
    dVector6& operator=(const dVector6& rhs)
    {
        memmove(data(), rhs.const_data(), 6*sizeof(double));
        return *this;
    };


//    double norm(void) const  /// Frobenius norm
//    {
//        double tmp = storage[0] * storage[0]
//                   + storage[1] * storage[1]
//                   + storage[2] * storage[2]
//                   + storage[3] * storage[3] * 2.0
//                   + storage[4] * storage[4] * 2.0
//                   + storage[5] * storage[5] * 2.0;
//        return sqrt(tmp);
//    };

    dVector6 rotated(const dMatrix3x3& RotationMatrix) const
    {
        double In[3][3];
        double Out[3][3];

        In[0][0] = storage[0];
        In[0][1] = storage[5];
        In[0][2] = storage[4];
        In[1][0] = storage[5];
        In[1][1] = storage[1];
        In[1][2] = storage[3];
        In[2][0] = storage[4];
        In[2][1] = storage[3];
        In[2][2] = storage[2];

        for(int p = 0; p < 3; ++p)
        for(int q = 0; q < 3; ++q)
        {
            Out[p][q] = 0;
            for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 3; ++j)
            {
                Out[p][q] += RotationMatrix(p,i)*In[i][j]*RotationMatrix(q,j);
            }
        }
        dVector6 OUT;

        OUT[0] = Out[0][0];
        OUT[5] = Out[0][1];
        OUT[4] = Out[0][2];
        OUT[1] = Out[1][1];
        OUT[3] = Out[1][2];
        OUT[2] = Out[2][2];

        return OUT;
    };
    dVector6& rotate(dMatrix3x3 RotationMatrix)
    {
        double In[3][3];
        double Out[3][3];

        In[0][0] = storage[0];
        In[0][1] = storage[5];
        In[0][2] = storage[4];
        In[1][0] = storage[5];
        In[1][1] = storage[1];
        In[1][2] = storage[3];
        In[2][0] = storage[4];
        In[2][1] = storage[3];
        In[2][2] = storage[2];

        for(int p = 0; p < 3; ++p)
        for(int q = 0; q < 3; ++q)
        {
            Out[p][q] = 0;
            for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 3; ++j)
            {
                Out[p][q] += RotationMatrix(p,i)*In[i][j]*RotationMatrix(q,j);
            }
        }

        storage[0] = Out[0][0];
        storage[5] = Out[0][1];
        storage[4] = Out[0][2];
        storage[1] = Out[1][1];
        storage[3] = Out[1][2];
        storage[2] = Out[2][2];

        return *this;
    };
    double* data(void)
    {
        return storage;
    };
    const double* const_data(void) const
    {
        return storage;
    };
    dVector6& H_product(dVector6 vector)
    {
        storage[0] *= vector[0];
        storage[1] *= vector[1];
        storage[2] *= vector[2];
        storage[3] *= vector[3];
        storage[4] *= vector[4];
        storage[5] *= vector[5];

        return *this;
    };

    dMatrix3x3 sym_tensor(void) const
    {
        dMatrix3x3 tmp;
        tmp(0,0) = storage[0];
        tmp(0,1) = storage[5];
        tmp(0,2) = storage[4];
        tmp(1,0) = storage[5];
        tmp(1,1) = storage[1];
        tmp(1,2) = storage[3];
        tmp(2,0) = storage[4];
        tmp(2,1) = storage[3];
        tmp(2,2) = storage[2];
        return tmp;
    };

    dMatrix3x3 skew_tensor(void) const
    {
        dMatrix3x3 tmp;
        tmp(0,0) =  storage[0];
        tmp(0,1) =  storage[5];
        tmp(0,2) =  storage[4];
        tmp(1,0) = -storage[5];
        tmp(1,1) =  storage[1];
        tmp(1,2) =  storage[3];
        tmp(2,0) = -storage[4];
        tmp(2,1) = -storage[3];
        tmp(2,2) =  storage[2];
        return tmp;
    };

    std::string print(void) const
    {
        std::stringstream out;
        out << "< | ";
        for(int i = 0; i < 6; i++)
        {
            out << storage[i]<< " " << " | ";
        }
        out << " >";
        return out.str();
    };

/*#ifdef _OPENMP
	static dVector6 maxCompdVector6(dVector6 vReduced, dVector6 vNew)
	{
		for (int n = 0; n < 6; n++)
			if (vNew[n]>vReduced[n])
			{
				vReduced[n] = vNew[n];
			}
		return vReduced;
	}
#pragma omp declare reduction(dVector6MaxComp:dVector6:omp_out=maxCompdVector6(omp_out,omp_in)) \
	  initializer(omp_priv=dVector6({-DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX, -DBL_MAX}))
#endif*/

 protected:
    double storage[6];
 private:
};

class dMatrix6x6
{
 public:
    double& operator()(const int i, const int j)
    {
#ifdef DEBUG
        if(i > 5 or j > 5)
        {
            std::cout << "Error in dMatrix6x6::operator()\n"
                 << "Access beyond storage range. (i, j) = "
                 << i <<", "<< j << " > (5, 5)"
                 << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i][j];
    };
    double const& operator()(const int i, const int j) const
    {
#ifdef DEBUG
        if(i > 5 or j > 5)
        {
            std::cout << "Error in dMatrix6x6::operator()\n"
                 << "Access beyond storage range. (i, j) = "
                 << i <<", "<< j << " > (5, 5)"
                 << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i][j];
    };
    void set_to_zero(void)
    {
        memset(storage, 0, 36*sizeof(double));
    };

    double norm(void) const
    {
        double tmp = 0.0;
        for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
        {
            tmp += storage[i][j]*storage[i][j];
        }
        return sqrt(tmp);
    };
    /*double det(void) const
    {
        double determinant = 0.0;

        for (int i = 0; i < 6; i++)
        {
            double line_product = 1.0;
            for (int j = 0; j < 6; j++)
            {
                line_product *= storage[(i+j)%6][j];
            }
            determinant += line_product;
        }

        for (int i = 0; i < 6; i++)
        {
            double line_product = 1.0;
            for (int j = 0; j < 6; j++)
            {
                line_product *= storage[(i-j+6)%6][j];
            }
            determinant -= line_product;
        }
        std::cout << determinant << std::endl;
        return determinant;
    };

    bool is_singular(void)
    {
        if(fabs(det()) > DBL_EPSILON)
        {
            return false;
        }
        else
        {
            return true;
        }
    };*/

    dMatrix6x6 operator*(const double m) const
    {
        dMatrix6x6 tmp;
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            tmp(i,j) = storage[i][j]*m;
        }
        return tmp;
    };
    dMatrix6x6& operator*=(const double m)
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            storage[i][j] *= m;
        }
        return *this;
    };
    /*dVector6 operator*(dVector6 rhs) const
    {
        dVector6 tmp;
        tmp.set_to_zero();
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            tmp[i] += storage[i][j]*rhs[j];
        }
        return tmp;
    };*/
    dMatrix6x6 operator*(const dMatrix6x6& rhs) const
    {
        dMatrix6x6 tmp;
        tmp.set_to_zero();
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        for(int k = 0; k < 6; k++)
        {
            tmp(i,j) += storage[i][k]*rhs(k,j);
        }
        return tmp;
    };

    vStrain operator*(const vStress& rhs) const;

    vStress operator*(const vStrain& rhs) const;

    dMatrix6x6 operator+(const dMatrix6x6& rhs) const
    {
        dMatrix6x6 tmp;
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            tmp(i,j) = storage[i][j] + rhs(i,j);
        }
        return tmp;
    };
    dMatrix6x6 operator-(const dMatrix6x6& rhs) const
    {
        dMatrix6x6 tmp;
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            tmp(i,j) = storage[i][j] - rhs(i,j);
        }
        return tmp;
    };
    dMatrix6x6& operator+=(const dMatrix6x6& rhs)
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            storage[i][j] += rhs(i,j);
        }
        return *this;
    };
    dMatrix6x6& operator/=(const double m)
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            storage[i][j] /= m;
        }
        return *this;
    };
    dMatrix6x6& operator-=(const dMatrix6x6& rhs)
    {
        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            storage[i][j] -= rhs(i,j);
        }
        return *this;
    };
    dMatrix6x6& operator=(const dMatrix6x6& rhs)
    {
        memmove(data(), rhs.const_data(), 36*sizeof(double));
        return *this;
    };
    dMatrix6x6& invert(void)
    {
        double Out[6][6];

        int indxc[6];
        int indxr[6];
        int ipiv[6] = {0, 0, 0, 0, 0, 0};
        int icol = 0;
        int irow = 0;
        double pivinv;
        double dum;
        double Uni[6][6] = {{1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
                            {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};


        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            Out[i][j] = storage[i][j];
        }

        for(int i = 0; i < 6; i++)
        {
            double big = 0.0;
            for(int j = 0; j < 6; j++)
            if(ipiv[j] != 1)
            for(int k = 0; k < 6; k++)
            {
                if(ipiv[k] == 0)
                {
                    if(fabs(Out[j][k]) >= big)
                    {
                        big = fabs(Out[j][k]);
                        irow = j;
                        icol = k;
                    };
                }
                else if (ipiv[k] > 1)
                {
                    std::cout << "dMatrix6x6: Can Not Compute Inverse Matrix."
                              << "Matrix:\n" << this->print()
                              << "is Singular 1!!!" << std::endl;
                    exit(1);
                }
            };
            ++(ipiv[icol]);
            if(irow != icol)
            {
                for (int l = 0; l < 6; l++)
                {
                    double temp = Out[irow][l];
                    Out[irow][l] = Out[icol][l];
                    Out[icol][l] = temp;
                };
                for (int l = 0; l < 6; l++)
                {
                    double temp = Uni[irow][l];
                    Uni[irow][l] = Uni[icol][l];
                    Uni[icol][l] = temp;
                };
            };
            indxr[i] = irow;
            indxc[i] = icol;
            if (fabs(Out[icol][icol]) <= DBL_EPSILON)
            {
                std::cout << "dMatrix6x6: Can Not Compute Inverse Matrix. Matrix:\n"<< this->print() << "is Singular 2!!!" << std::endl;
                exit(2);
            }
            pivinv = 1.0/Out[icol][icol];
            Out[icol][icol] = 1.0;
            for(int l = 0; l < 6; l++) Out[icol][l] *= pivinv;
            for(int l = 0; l < 6; l++) Uni[icol][l] *= pivinv;
            for(int ll = 0; ll < 6; ll++)
            if(ll != icol)
            {
                dum = Out[ll][icol];
                Out[ll][icol] = 0.0;
                for(int l = 0; l < 6; l++) Out[ll][l] -= Out[icol][l]*dum;
                for(int l = 0; l < 6; l++) Uni[ll][l] -= Uni[icol][l]*dum;
            }
        }
        for(int l = 5; l >= 0; l--)
        {
            if(indxr[l] != indxc[l])
            for(int k = 0; k < 6; k++)
            {
                double temp = Out[k][indxr[l]];
                Out[k][indxr[l]] = Out[k][indxc[l]];
                Out[k][indxc[l]] = temp;
            };
        }

        memmove(storage, Out, 36*sizeof(double));
        return *this;
    };

    dMatrix6x6 inverted(void) const
    {
        double Out[6][6];

        int indxc[6];
        int indxr[6];
        int ipiv[6] = {0, 0, 0, 0, 0, 0};
        int icol = 0;
        int irow = 0;
        double pivinv;
        double dum;
        double Uni[6][6] = {{1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                            {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                            {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
                            {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};

        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            Out[i][j] = storage[i][j];
        }

        for(int i = 0; i < 6; i++)
        {
            double big = 0.0;
            for(int j = 0; j < 6; j++)
            if(ipiv[j] != 1)
            for(int k = 0; k < 6; k++)
            {
                if(ipiv[k] == 0)
                {
                    if(fabs(Out[j][k]) >= big)
                    {
                        big = fabs(Out[j][k]);
                        irow = j;
                        icol = k;
                    };
                }
                else if (ipiv[k] > 1)
                {
                    std::cout << "dMatrix6x6: Can Not Compute Inverse Matrix.\n"
                              << this->print() << "Matrix is Singular 2!!!"
                              << std::endl;
                    exit(1);
                }
            };
            ++(ipiv[icol]);
            if(irow != icol)
            {
                for (int l = 0; l < 6; l++)
                {
                    double temp = Out[irow][l];
                    Out[irow][l] = Out[icol][l];
                    Out[icol][l] = temp;
                };
                for (int l = 0; l < 6; l++)
                {
                    double temp = Uni[irow][l];
                    Uni[irow][l] = Uni[icol][l];
                    Uni[icol][l] = temp;
                };
            };
            indxr[i] = irow;
            indxc[i] = icol;
            if (fabs(Out[icol][icol]) <= DBL_EPSILON)
            {
                std::cout << "dMatrix6x6:\n" << this->print()
                          << " Is not Invertible!" << std::endl;
                exit(2);
            }
            pivinv = 1.0/Out[icol][icol];
            Out[icol][icol] = 1.0;
            for(int l = 0; l < 6; l++) Out[icol][l] *= pivinv;
            for(int l = 0; l < 6; l++) Uni[icol][l] *= pivinv;
            for(int ll = 0; ll < 6; ll++)
            if(ll != icol)
            {
                dum = Out[ll][icol];
                Out[ll][icol] = 0.0;
                for(int l = 0; l < 6; l++) Out[ll][l] -= Out[icol][l]*dum;
                for(int l = 0; l < 6; l++) Uni[ll][l] -= Uni[icol][l]*dum;
            }
        }
        for(int l = 5; l >= 0; l--)
        {
            if(indxr[l] != indxc[l])
            for(int k = 0; k < 6; k++)
            {
                double temp = Out[k][indxr[l]];
                Out[k][indxr[l]] = Out[k][indxc[l]];
                Out[k][indxc[l]] = temp;
            };
        }

        dMatrix6x6 tmp;
        tmp.set_to_zero();
        memmove(tmp.data(), Out, 36*sizeof(double));
        return tmp;
    };

    dMatrix6x6 set_to_unity(void)
    {
        double Uni[6][6] = {{1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                             {0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
                             {0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
                             {0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
                             {0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
                             {0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};

        dMatrix6x6 tmp;
        memmove(tmp.data(), Uni, 36*sizeof(double));
        return tmp;
    };

    dMatrix6x6& transpose(void)
    {
        double tmp[6][6];
        memset(tmp, 0, 36*sizeof(double));

        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            tmp[i][j] = storage[j][i];
        }
        memmove(storage, tmp, 36*sizeof(double));
        return *this;
    };

    dMatrix6x6 transposed(void) const
    {
        dMatrix6x6 tmp;

        for(int i = 0; i < 6; i++)
        for(int j = 0; j < 6; j++)
        {
            tmp(i,j) = storage[j][i];
        }
        return tmp;
    };

    dMatrix6x6& rotate(const dMatrix3x3& RotationMatrix)
    {
        double In[3][3][3][3];
        double Out[3][3][3][3];

        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        for(int k = 0; k < 3; k++)
        for(int l = 0; l < 3; l++)
        {
            In[i][j][k][l] = tensor(i,j,k,l);
        }

        for(int m = 0; m < 3; m++)
        for(int n = 0; n < 3; n++)
        for(int p = 0; p < 3; p++)
        for(int q = 0; q < 3; q++)
        {
            Out[m][n][p][q] = 0.0;

            for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
            for(int l = 0; l < 3; l++)
            {
//                Former version
//                Out[m][n][p][q] += RotationMatrix(i,m)*
//                                   RotationMatrix(j,n)*
//                                   In[i][j][k][l]*
//                                   RotationMatrix(k,p)*
//                                   RotationMatrix(l,q);
//              New version (active rotation)
                Out[m][n][p][q] += RotationMatrix(m,i)*
                                   RotationMatrix(n,j)*
                                   In[i][j][k][l]*
                                   RotationMatrix(p,k)*
                                   RotationMatrix(q,l);
            }
        }
        int VoigtIndex[6][2] = {{0,0},{1,1},{2,2},{1,2},{0,2},{0,1}};
        for(int m = 0; m < 6; m++)
        for(int n = 0; n < 6; n++)
        {
            int i = VoigtIndex[m][0];
            int j = VoigtIndex[m][1];
            int k = VoigtIndex[n][0];
            int l = VoigtIndex[n][1];

            storage[m][n] = Out[i][j][k][l];
        }
        return *this;
    };
    dMatrix6x6 rotated(const dMatrix3x3& RotationMatrix)
    {
        double In[3][3][3][3];
        double Out[3][3][3][3];

        for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++)
        for(int k = 0; k < 3; k++)
        for(int l = 0; l < 3; l++)
        {
            In[i][j][k][l] = tensor(i,j,k,l);
        }

        for(int m = 0; m < 3; m++)
        for(int n = 0; n < 3; n++)
        for(int p = 0; p < 3; p++)
        for(int q = 0; q < 3; q++)
        {
            Out[m][n][p][q] = 0.0;

            for(int i = 0; i < 3; i++)
            for(int j = 0; j < 3; j++)
            for(int k = 0; k < 3; k++)
            for(int l = 0; l < 3; l++)
            {
//                Former version
//                Out[m][n][p][q] += RotationMatrix(i,m)*
//                                   RotationMatrix(j,n)*
//                                   In[i][j][k][l]*
//                                   RotationMatrix(k,p)*
//                                   RotationMatrix(l,q);
//              New version (active rotation)
                Out[m][n][p][q] += RotationMatrix(m,i)*
                                   RotationMatrix(n,j)*
                                   In[i][j][k][l]*
                                   RotationMatrix(p,k)*
                                   RotationMatrix(q,l);
            }
        }
        dMatrix6x6 OUT;
        int VoigtIndex[6][2] = {{0,0},{1,1},{2,2},{1,2},{0,2},{0,1}};

        for(int m = 0; m < 6; m++)
        for(int n = 0; n < 6; n++)
        {
            int i = VoigtIndex[m][0];
            int j = VoigtIndex[m][1];
            int k = VoigtIndex[n][0];
            int l = VoigtIndex[n][1];

            OUT(m,n) = Out[i][j][k][l];
        }
        return OUT;
    };

    std::string print(void) const
    {
        std::stringstream out;
        for(int i = 0; i < 6; i++)
        {
            out << "||" << std::setprecision(6)
                        << std::setw(8) << storage[i][0] << " "
                        << std::setw(8) << storage[i][1] << " "
                        << std::setw(8) << storage[i][2] << " "
                        << std::setw(8) << storage[i][3] << " "
                        << std::setw(8) << storage[i][4] << " "
                        << std::setw(8) << storage[i][5] << "||\n";
        }
        return out.str();
    };
    const double& tensor(const int i, const int j, const int k, const int l) const
    {
        return storage[(i==j)?(i):(6-(i+j))][(k==l)?(k):(6-(k+l))];
    };
    double* data(void)
    {
        return &storage[0][0];
    };
    const double* const_data(void) const
    {
        return &storage[0][0];
    };
 protected:
 private:

    double storage[6][6];
};

inline dMatrix6x6 dMatrix3x3::outer(const dMatrix3x3& rhs) const
{
    dMatrix6x6 tmp;
    tmp.set_to_zero();

    tmp(0,0) += storage[0][0]*rhs(0, 0);
    tmp(0,5) += storage[0][0]*rhs(0, 1);
    tmp(0,4) += storage[0][0]*rhs(0, 2);
    tmp(0,5) += storage[0][0]*rhs(1, 0);
    tmp(0,1) += storage[0][0]*rhs(1, 1);
    tmp(0,3) += storage[0][0]*rhs(1, 2);
    tmp(0,4) += storage[0][0]*rhs(2, 0);
    tmp(0,3) += storage[0][0]*rhs(2, 1);
    tmp(0,2) += storage[0][0]*rhs(2, 2);
    tmp(5,0) += storage[0][1]*rhs(0, 0);
    tmp(5,5) += storage[0][1]*rhs(0, 1);
    tmp(5,4) += storage[0][1]*rhs(0, 2);
    tmp(5,5) += storage[0][1]*rhs(1, 0);
    tmp(5,1) += storage[0][1]*rhs(1, 1);
    tmp(5,3) += storage[0][1]*rhs(1, 2);
    tmp(5,4) += storage[0][1]*rhs(2, 0);
    tmp(5,3) += storage[0][1]*rhs(2, 1);
    tmp(5,2) += storage[0][1]*rhs(2, 2);
    tmp(4,0) += storage[0][2]*rhs(0, 0);
    tmp(4,5) += storage[0][2]*rhs(0, 1);
    tmp(4,4) += storage[0][2]*rhs(0, 2);
    tmp(4,5) += storage[0][2]*rhs(1, 0);
    tmp(4,1) += storage[0][2]*rhs(1, 1);
    tmp(4,3) += storage[0][2]*rhs(1, 2);
    tmp(4,4) += storage[0][2]*rhs(2, 0);
    tmp(4,3) += storage[0][2]*rhs(2, 1);
    tmp(4,2) += storage[0][2]*rhs(2, 2);
    tmp(5,0) += storage[1][0]*rhs(0, 0);
    tmp(5,5) += storage[1][0]*rhs(0, 1);
    tmp(5,4) += storage[1][0]*rhs(0, 2);
    tmp(5,5) += storage[1][0]*rhs(1, 0);
    tmp(5,1) += storage[1][0]*rhs(1, 1);
    tmp(5,3) += storage[1][0]*rhs(1, 2);
    tmp(5,4) += storage[1][0]*rhs(2, 0);
    tmp(5,3) += storage[1][0]*rhs(2, 1);
    tmp(5,2) += storage[1][0]*rhs(2, 2);
    tmp(1,0) += storage[1][1]*rhs(0, 0);
    tmp(1,5) += storage[1][1]*rhs(0, 1);
    tmp(1,4) += storage[1][1]*rhs(0, 2);
    tmp(1,5) += storage[1][1]*rhs(1, 0);
    tmp(1,1) += storage[1][1]*rhs(1, 1);
    tmp(1,3) += storage[1][1]*rhs(1, 2);
    tmp(1,4) += storage[1][1]*rhs(2, 0);
    tmp(1,3) += storage[1][1]*rhs(2, 1);
    tmp(1,2) += storage[1][1]*rhs(2, 2);
    tmp(3,0) += storage[1][2]*rhs(0, 0);
    tmp(3,5) += storage[1][2]*rhs(0, 1);
    tmp(3,4) += storage[1][2]*rhs(0, 2);
    tmp(3,5) += storage[1][2]*rhs(1, 0);
    tmp(3,1) += storage[1][2]*rhs(1, 1);
    tmp(3,3) += storage[1][2]*rhs(1, 2);
    tmp(3,4) += storage[1][2]*rhs(2, 0);
    tmp(3,3) += storage[1][2]*rhs(2, 1);
    tmp(3,2) += storage[1][2]*rhs(2, 2);
    tmp(4,0) += storage[2][0]*rhs(0, 0);
    tmp(4,5) += storage[2][0]*rhs(0, 1);
    tmp(4,4) += storage[2][0]*rhs(0, 2);
    tmp(4,5) += storage[2][0]*rhs(1, 0);
    tmp(4,1) += storage[2][0]*rhs(1, 1);
    tmp(4,3) += storage[2][0]*rhs(1, 2);
    tmp(4,4) += storage[2][0]*rhs(2, 0);
    tmp(4,3) += storage[2][0]*rhs(2, 1);
    tmp(4,2) += storage[2][0]*rhs(2, 2);
    tmp(3,0) += storage[2][1]*rhs(0, 0);
    tmp(3,5) += storage[2][1]*rhs(0, 1);
    tmp(3,4) += storage[2][1]*rhs(0, 2);
    tmp(3,5) += storage[2][1]*rhs(1, 0);
    tmp(3,1) += storage[2][1]*rhs(1, 1);
    tmp(3,3) += storage[2][1]*rhs(1, 2);
    tmp(3,4) += storage[2][1]*rhs(2, 0);
    tmp(3,3) += storage[2][1]*rhs(2, 1);
    tmp(3,2) += storage[2][1]*rhs(2, 2);
    tmp(2,0) += storage[2][2]*rhs(0, 0);
    tmp(2,5) += storage[2][2]*rhs(0, 1);
    tmp(2,4) += storage[2][2]*rhs(0, 2);
    tmp(2,5) += storage[2][2]*rhs(1, 0);
    tmp(2,1) += storage[2][2]*rhs(1, 1);
    tmp(2,3) += storage[2][2]*rhs(1, 2);
    tmp(2,4) += storage[2][2]*rhs(2, 0);
    tmp(2,3) += storage[2][2]*rhs(2, 1),
    tmp(2,2) += storage[2][2]*rhs(2, 2);

    return tmp;
}

inline dVector6 dMatrix3x3::VoigtVector() const
{
    dVector6 tmp;
    tmp[0] = storage[0][0];
    tmp[1] = storage[1][1];
    tmp[2] = storage[2][2];
    tmp[3] = storage[1][2];
    tmp[4] = storage[0][2];
    tmp[5] = storage[0][1];
    return tmp;
}

class vStress: public dVector6
{
    /*
     * Stress Voigt vector
     */
public:
    vStress()
    {
        memset(data(), 0, 6*sizeof(double));
    };
    /*vStress(vStress& rhs)
    {
        memmove(data(), rhs.data(), 6*sizeof(double));
    };*/
    vStress(const dVector6& rhs)
    {
        memmove(data(), rhs.const_data(), 6*sizeof(double));
    };
    vStress& operator=(const dVector6& rhs)
    {
        memmove(data(), rhs.const_data(), 6*sizeof(double));
        return *this;
    };
    /*vStress& operator=(vStress& rhs)
    {
        memmove(data(), rhs.data(), 6*sizeof(double));
        return *this;
    };*/
    vStress operator-(const vStress& rhs) const
    {
        vStress tmp;
        tmp[0] = storage[0] - rhs[0];
        tmp[1] = storage[1] - rhs[1];
        tmp[2] = storage[2] - rhs[2];
        tmp[3] = storage[3] - rhs[3];
        tmp[4] = storage[4] - rhs[4];
        tmp[5] = storage[5] - rhs[5];
        return tmp;
    };
    double norm(void) const  /// Frobenius norm
    {
        double tmp = storage[0] * storage[0]
                   + storage[1] * storage[1]
                   + storage[2] * storage[2]
                   + storage[3] * storage[3] * 2.0
                   + storage[4] * storage[4] * 2.0
                   + storage[5] * storage[5] * 2.0;
        return sqrt(tmp);
    };

    double doublecontract(const vStress& Bstress) const  // "double-dot product"
    {
        double tmp =
        storage[0] * Bstress[0] +
        storage[1] * Bstress[1] +
        storage[2] * Bstress[2] +
        2.0*storage[3] * Bstress[3] +
        2.0*storage[4] * Bstress[4] +
        2.0*storage[5] * Bstress[5];
        return tmp;
    };

    double doublecontract(const dVector6& symTensorV) const  // "double-dot product"
    {
        double tmp =
        storage[0] * symTensorV[0] +
        storage[1] * symTensorV[1] +
        storage[2] * symTensorV[2] +
        2.0*storage[3] * symTensorV[3] +
        2.0*storage[4] * symTensorV[4] +
        2.0*storage[5] * symTensorV[5];
        return tmp;
    };

    inline double Pressure() const
    {
        return -(storage[0] + storage[1] + storage[2])/3.0;
    };

    inline double Trace() const
    {
        return storage[0] + storage[1] + storage[2];
    };

    inline double Determinant() const
    {
        return storage[0]*storage[1]*storage[2]-
            storage[0]*storage[3]*storage[3] -
            storage[1]*storage[4]*storage[4] -
            storage[2]*storage[5]*storage[5] +
            storage[3]*storage[4]*storage[5]*2;
    };

    dVector3 Invariants(void) const
    {
        const double I1 = Trace();
        const double I2 = 0.5*(Trace()*Trace() -
                (storage[0]*storage[0] +
                 storage[1]*storage[1] +
                 storage[2]*storage[2] +
                 storage[3]*storage[3]*2 +
                 storage[4]*storage[4]*2 +
                 storage[5]*storage[5]*2));
        const double I3 = Determinant();
        return dVector3({I1,I2,I3});
    }

    double Mises(void) const
    {
        double vMises = 0.0;
        vMises = (storage[0]-storage[1])*(storage[0]-storage[1]) +
                 (storage[1]-storage[2])*(storage[1]-storage[2]) +
                 (storage[2]-storage[0])*(storage[2]-storage[0]) +
                 6.*(storage[3]*storage[3]+storage[4]*storage[4]+storage[5]*storage[5]);

        return sqrt(0.5*vMises);
    }

    vStress rotated(const dMatrix3x3& RotationMatrix) const
    {
        double In[3][3];
        double Out[3][3];

        In[0][0] = storage[0];
        In[0][1] = storage[5];
        In[0][2] = storage[4];
        In[1][0] = storage[5];
        In[1][1] = storage[1];
        In[1][2] = storage[3];
        In[2][0] = storage[4];
        In[2][1] = storage[3];
        In[2][2] = storage[2];

        for(int p = 0; p < 3; ++p)
        for(int q = 0; q < 3; ++q)
        {
            Out[p][q] = 0;
            for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 3; ++j)
            {
                Out[p][q] += RotationMatrix(p,i)*In[i][j]*RotationMatrix(q,j);
            }
        }
        vStress OUT;

        OUT[0] = Out[0][0];
        OUT[5] = Out[0][1];
        OUT[4] = Out[0][2];
        OUT[1] = Out[1][1];
        OUT[3] = Out[1][2];
        OUT[2] = Out[2][2];

        return OUT;
    };
    vStress& rotate(const dMatrix3x3& RotationMatrix)
    {
        double In[3][3];
        double Out[3][3];

        In[0][0] = storage[0];
        In[0][1] = storage[5];
        In[0][2] = storage[4];
        In[1][0] = storage[5];
        In[1][1] = storage[1];
        In[1][2] = storage[3];
        In[2][0] = storage[4];
        In[2][1] = storage[3];
        In[2][2] = storage[2];

        for(int p = 0; p < 3; ++p)
        for(int q = 0; q < 3; ++q)
        {
            Out[p][q] = 0;
            for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 3; ++j)
            {
                Out[p][q] += RotationMatrix(p,i)*In[i][j]*RotationMatrix(q,j);
            }
        }

        storage[0] = Out[0][0];
        storage[5] = Out[0][1];
        storage[4] = Out[0][2];
        storage[1] = Out[1][1];
        storage[3] = Out[1][2];
        storage[2] = Out[2][2];

        return *this;
    };

    double get_tensor(const int i, const int j) const
    {
        return storage[(i==j)?(i):(6-(i+j))];
    };

    dMatrix3x3 tensor(void) const
    {
        dMatrix3x3 tmp;
        tmp(0,0) = storage[0];
        tmp(0,1) = storage[5];
        tmp(0,2) = storage[4];
        tmp(1,0) = storage[5];
        tmp(1,1) = storage[1];
        tmp(1,2) = storage[3];
        tmp(2,0) = storage[4];
        tmp(2,1) = storage[3];
        tmp(2,2) = storage[2];
        return tmp;
    };
protected:
private:
};

class vStrain: public dVector6
{
    /*
     *  This is a special version of the six component Voigt vector:
     *  the off diagonal elements of the strain matrix are multiplied by 2(!)
     *  while transforming the 3x3 strain matrix to a 6-component vector,
     *  and have to be divided by two if they are transformed to a 3x3 matrix
     *  again. Therefore the transformation functions and some other functions
     *  have to be specified in different way then for the stress-type
     *  Voigt vector. The differences are specified in the functions.
     */
 public:
    vStrain()
    {
        memset(data(), 0, 6*sizeof(double));
    };
    vStrain(const dVector6& rhs)
    {
        memmove(data(), rhs.const_data(), 6*sizeof(double));
    };
    vStrain& operator=(const dVector6& rhs)
    {
        memmove(data(), rhs.const_data(), 6*sizeof(double));
        return *this;
    };
    bool operator==(const dVector6& rhs)
    {
        for (int i = 0; i < 6; i++)
        {
            if (storage[i] != rhs[i]){return false;};
        }
        return true;
    };
    vStrain operator-(const vStrain& rhs) const
    {
        vStrain tmp;
        tmp[0] = storage[0] - rhs[0];
        tmp[1] = storage[1] - rhs[1];
        tmp[2] = storage[2] - rhs[2];
        tmp[3] = storage[3] - rhs[3];
        tmp[4] = storage[4] - rhs[4];
        tmp[5] = storage[5] - rhs[5];
        return tmp;
    };
    double norm(void) const  /// Frobenius norm
    {
        /*
         * multiplication by 0.5 of the off diagonal elements before taking
         * square and doubling the result due to double appearance of the off
         * diagonal elements => 0.5^2 / 2 = 0.5
         */
        double tmp = storage[0] * storage[0]
                   + storage[1] * storage[1]
                   + storage[2] * storage[2]
                   + storage[3] * storage[3] * 0.5
                   + storage[4] * storage[4] * 0.5
                   + storage[5] * storage[5] * 0.5;
        return sqrt(tmp);
    };

    double doublecontract(const vStrain& Bstrain) const // "double-dot product"
    {
        double tmp =
        storage[0] * Bstrain[0] +
        storage[1] * Bstrain[1] +
        storage[2] * Bstrain[2] +
        2.0*(0.5*storage[3]) * (0.5*Bstrain[3]) +
        2.0*(0.5*storage[4]) * (0.5*Bstrain[4]) +
        2.0*(0.5*storage[5]) * (0.5*Bstrain[5]);
        return tmp;
    };

    vStrain rotated(const dMatrix3x3& RotationMatrix)
    {
        /*
         * multiplication by 0.5 of the off diagonal elements during translating
         * to a 3x3 matrix and multiplying by two during translating back to
         * a Voigt strain vector
         */
        double In[3][3];
        double Out[3][3];

        In[0][0] = storage[0];
        In[0][1] = storage[5]*0.5;
        In[0][2] = storage[4]*0.5;
        In[1][0] = storage[5]*0.5;
        In[1][1] = storage[1];
        In[1][2] = storage[3]*0.5;
        In[2][0] = storage[4]*0.5;
        In[2][1] = storage[3]*0.5;
        In[2][2] = storage[2];

        for(int p = 0; p < 3; ++p)
        for(int q = 0; q < 3; ++q)
        {
            Out[p][q] = 0;
            for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 3; ++j)
            {
                Out[p][q] += RotationMatrix(p,i)*In[i][j]*RotationMatrix(q,j);
            }
        }
        vStrain OUT;

        OUT[0] = Out[0][0];
        OUT[5] = Out[0][1]*2.0;
        OUT[4] = Out[0][2]*2.0;
        OUT[1] = Out[1][1];
        OUT[3] = Out[1][2]*2.0;
        OUT[2] = Out[2][2];

        return OUT;
    };

    vStrain& rotate(const dMatrix3x3& RotationMatrix)
    {
        /*
         * multiplication by 0.5 of the off diagonal elements during translating
         * to a 3x3 matrix and multiplying by two during translating back to
         * a Voigt strain vector
         */
        double In[3][3];
        double Out[3][3];

        In[0][0] = storage[0];
        In[0][1] = storage[5]*0.5;
        In[0][2] = storage[4]*0.5;
        In[1][0] = storage[5]*0.5;
        In[1][1] = storage[1];
        In[1][2] = storage[3]*0.5;
        In[2][0] = storage[4]*0.5;
        In[2][1] = storage[3]*0.5;
        In[2][2] = storage[2];

        for(int p = 0; p < 3; ++p)
        for(int q = 0; q < 3; ++q)
        {
            Out[p][q] = 0;
            for(int i = 0; i < 3; ++i)
            for(int j = 0; j < 3; ++j)
            {
                Out[p][q] += RotationMatrix(p,i)*In[i][j]*RotationMatrix(q,j);
            }
        }

        storage[0] = Out[0][0];
        storage[5] = Out[0][1]*2.0;
        storage[4] = Out[0][2]*2.0;
        storage[1] = Out[1][1];
        storage[3] = Out[1][2]*2.0;
        storage[2] = Out[2][2];

        return *this;
    };

    vStrain Ln() const;                                                     // Converts the strain to the logarithmic strain

    double get_tensor(const int i, const int j) const
    {
        /*
         * multiplication by 0.5 of the off diagonal elements during translating
         * to a 3x3 matrix
         */
        return (storage[(i==j)?(i):(6-(i+j))] * ((i==j)?(1.0):(0.5)));
    };

    dMatrix3x3 tensor(void) const
    {
        /*
         * multiplication by 0.5 of the off diagonal elements during translating
         * to a 3x3 matrix
         */
        dMatrix3x3 tmp;
        tmp(0,0) = storage[0];
        tmp(0,1) = storage[5]*0.5;
        tmp(0,2) = storage[4]*0.5;
        tmp(1,0) = storage[5]*0.5;
        tmp(1,1) = storage[1];
        tmp(1,2) = storage[3]*0.5;
        tmp(2,0) = storage[4]*0.5;
        tmp(2,1) = storage[3]*0.5;
        tmp(2,2) = storage[2];
        return tmp;
    };
};

inline vStrain dMatrix3x3::VoigtStrain() const
{
    vStrain tmp;
    tmp[0] = storage[0][0];
    tmp[1] = storage[1][1];
    tmp[2] = storage[2][2];
    tmp[3] = storage[1][2]*2.0;
    tmp[4] = storage[0][2]*2.0;
    tmp[5] = storage[0][1]*2.0;
    return tmp;
}

inline vStress dMatrix3x3::VoigtStress() const
{
    vStress tmp;
    tmp[0] = storage[0][0];
    tmp[1] = storage[1][1];
    tmp[2] = storage[2][2];
    tmp[3] = storage[1][2];
    tmp[4] = storage[0][2];
    tmp[5] = storage[0][1];
    return tmp;
}

inline vStrain dMatrix6x6::operator*(const vStress& rhs) const
{
    vStrain tmp;
    tmp.set_to_zero();
    for(int i = 0; i < 6; i++)
    for(int j = 0; j < 6; j++)
    {
        tmp[i] += storage[i][j]*rhs[j];
    }
    return tmp;
}

inline vStress dMatrix6x6::operator*(const vStrain& rhs) const
{
    vStress tmp;
    tmp.set_to_zero();
    for(int i = 0; i < 6; i++)
    for(int j = 0; j < 6; j++)
    {
        tmp[i] += storage[i][j]*rhs[j];
    }
    return tmp;
}

}// namespace opensim
#endif
