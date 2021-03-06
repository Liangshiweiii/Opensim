#ifndef QUATERNION_H
#define QUATERNION_H

#include "Tools/Includes.h"
namespace opensim
{

class Quaternion                                                                /// Stores and manages Quaternions
{
 public:

    Quaternion()                                                               /// Standard constructor
    {
        s = 0.0;
        v.set_to_zero();
        RotationMatrix.set_to_unity();
    }
    Quaternion(const Quaternion& rhS)
    {
        s = rhS.s;
        v = rhS.v;
        setRotationMatrix();
    }
    Quaternion(std::initializer_list<double> vecinit)
    {
#ifdef DEBUG
        if (vecinit.size() != 4)
        {
            std::cout << "Error in Quaternion::constructor()\n"
                      << "Initialization list size beyond storage range."
                      << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif

        int ii = 0;
        s = *vecinit.begin();
        for (auto it = vecinit.begin()+1; it != vecinit.end(); it++)
        {
            v[ii] = *it;
            ii += 1;
        }
        setRotationMatrix();
    }

    double& operator[](const int index)
    {
    #ifdef DEBUG
        if(index > 3 or index < 0)
        {
            std::cout << "Error in Quaternion::operator[]\n"
                      << "Access beyond storage range. i = "
                      << index << " > 3"
                      << "\nTerminating!!!" << std::endl;
            exit(13);
        }
    #endif
        if (index == 0){return s;};
        return v[index-1];
    }

    double const& operator[](const int index) const
    {
    #ifdef DEBUG
        if(index > 3 or index < 0)
        {
            std::cout << "Error in Quaternion::operator[]\n"
                      << "Access beyond storage range. i = "
                      << index << " > 3"
                      << "\nTerminating!!!" << std::endl;
            exit(13);
        }
    #endif
        if (index == 0){return s;};
        return v[index-1];
    }
    void set(const double sIn, const double x, const double y, const double z)
    {
        s = sIn;
        v[0] = x;
        v[1] = y;
        v[2] = z;
        setRotationMatrix();
    }
    
    void set(const double sIn, const dVector3 vIn)
    {
        s = sIn;
        v = vIn;
        setRotationMatrix();
    }
    
    void set(const dMatrix3x3& RotMatrix)
    {
#ifdef DEBUG
        if (RotMatrix.determinant() < 1.0 - 3.0*DBL_EPSILON or 
            RotMatrix.determinant() > 1.0 + 3.0*DBL_EPSILON)
        {
            std::cout << "Supplied matrix is not a proper rotation matrix.\n"
                      << "Orientation::operator=.\n"
                      << "Terminating!" << std::endl;
            exit(1);
        }
#endif
        // Found at http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
        RotationMatrix = RotMatrix;
        
        double trace = RotMatrix.trace();

        if(1.0 + trace > 0)
        {
            double ss = sqrt(trace + 1.0)*2.0;
            set(0.25*ss,
               (RotMatrix(2,1) - RotMatrix(1,2))/ss,
               (RotMatrix(0,2) - RotMatrix(2,0))/ss,
               (RotMatrix(1,0) - RotMatrix(0,1))/ss);
        }
        else
        {
            if (RotMatrix(0,0) > RotMatrix(1,1) && RotMatrix(0,0) > RotMatrix(2,2))
            {
                double ss = 2.0 * sqrt( 1.0 + RotMatrix(0,0) - RotMatrix(1,1) - RotMatrix(2,2));
                set((RotMatrix(2,1) - RotMatrix(1,2))/ss,
                     0.25 * ss,
                    (RotMatrix(0,1) + RotMatrix(1,0))/ss,
                    (RotMatrix(0,2) + RotMatrix(2,0))/ss);
            }
            else if (RotMatrix(1,1) > RotMatrix(2,2))
            {
                double ss = 2.0 * sqrt(1.0 + RotMatrix(1,1) - RotMatrix(0,0) - RotMatrix(2,2));
                set((RotMatrix(0,2) - RotMatrix(2,0))/ss,
                    (RotMatrix(0,1) + RotMatrix(1,0))/ss,
                     0.25 * ss,
                    (RotMatrix(1,2) + RotMatrix(2,1))/ss);
            }
            else
            {
                double ss = 2.0 * sqrt(1.0 + RotMatrix(2,2) - RotMatrix(0,0) - RotMatrix(1,1));
                set((RotMatrix(1,0) - RotMatrix(0,1))/ss,
                    (RotMatrix(0,2) + RotMatrix(2,0))/ss,
                    (RotMatrix(1,2) + RotMatrix(2,1))/ss,
                     0.25 * ss);
            }
        }
    }
    
    void set_entry(const int idx, const double val)
    {
        if (idx == 0)
        {
            s = val;
        }
        else
        {
            v[idx] = val;
        }
        setRotationMatrix();
    }
    void set_to_zero()
    {
        s    = 0.0;
        v[0] = 0.0;
        v[1] = 0.0;
        v[2] = 0.0;
        RotationMatrix.set_to_unity();
    }
    Quaternion& operator=(const Quaternion rhS)
    {
        s = rhS.s;
        v = rhS.v;
        setRotationMatrix();
        return *this;
    }
    Quaternion  operator+(const Quaternion rhS) const
    {
        Quaternion result;

        result.s = s + rhS.s;
        result.v = v + rhS.v;
        //result.setRotationMatrix();
        return result;
    }
    Quaternion& operator+=(const Quaternion rhS)
    {
        s = s + rhS.s;
        v = v + rhS.v;
        setRotationMatrix();
        return *this;
    }
    Quaternion  operator-(const Quaternion rhS) const
    {
        Quaternion result;
        result.s = s - rhS.s;
        result.v = v - rhS.v;
        //result.setRotationMatrix();
        return result;
    }
    Quaternion& operator-=(const Quaternion rhS)
    {
        s = s - rhS.s;
        v = v - rhS.v;
        setRotationMatrix();
        return *this;
    }
    Quaternion  operator*(const double scalar) const
    {
        Quaternion result;
        result.s = s * scalar;
        result.v = v * scalar;
        //result.setRotationMatrix();
        return result;
    }
    Quaternion& operator*=(const double scalar)
    {
        s *= scalar;
        v *= scalar;
        setRotationMatrix();
        return *this;
    }
    Quaternion  operator*(const Quaternion rhS) const
    {
        Quaternion result;
        //result.set(s*rhS.s - v*rhS.v,v.cross(rhS.v) + rhS.v*s + v*rhS.s);

        result.s    =        s*rhS.s - v[1]*rhS.v[1] - v[2]*rhS.v[2] - v[3]*rhS.v[3];
        result.v[0] = rhS.v[1]*rhS.s +    s*rhS.v[1] - v[3]*rhS.v[2] + v[2]*rhS.v[3];
        result.v[1] = rhS.v[2]*rhS.s + v[3]*rhS.v[1] +   s *rhS.v[2] - v[1]*rhS.v[3];
        result.v[2] = rhS.v[3]*rhS.s - v[2]*rhS.v[1] + v[1]*rhS.v[2] +    s*rhS.v[3];
        //result.setRotationMatrix();
        return result;
    }
    Quaternion& operator*=(const Quaternion rhS)
    {
        s = s*rhS.s - v*rhS.v;
        v = v.cross(rhS.v) + rhS.v*s + v*rhS.s;
        setRotationMatrix();
        return *this;
    }
    Quaternion  operator/(const double divisor)
    {
        Quaternion result;
        result.s = s/divisor;
        result.v = v/divisor;
        //result.setRotationMatrix();
        return result;
    }
    Quaternion& operator/=(const double divisor)
    {
        s = s/divisor;
        v = v/divisor;
        setRotationMatrix();
        return *this;
    }
    double length() const
    {
        return sqrt(s*s + v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }
    Quaternion& normalize()
    {
        *this /= sqrt(s*s + v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        return *this;
    }
    Quaternion  normalized() const
    {
        Quaternion result;
        double norm = length();
        result.s = s/norm;
        result.v = v/norm;
        return result;
    }
    Quaternion& conjugate()
    {
        v = v*(-1.0);
        setRotationMatrix();
        return *this;
    }
    Quaternion  conjugated() const
    {
        Quaternion result;
        result.s = s;
        result.v = v*(-1.0);
        return result;
    }
    Quaternion& invert()
    {
        // if |q|=1 then q inverse = q conjugate
        v = v*(-1.0);
        *this /= (s*s + v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        return *this;
    }
    Quaternion  inverted() const
    {
        Quaternion result = *this;
        return result = result.conjugate()/(s*s+v*v);
    }

    void setRotationMatrix(void);
    
    //EulerAngles getEulerAngles(const EulerConvention EConvention) const;

    static Quaternion lerp(const Quaternion& rhSQ1, const Quaternion& rhSQ2,
                                                            const double t);
    static Quaternion slerp(const Quaternion& a, const Quaternion& b,
                                                            const double t);
    std::string print(void) const;
    std::string print_matrix(void) const;
    
    dMatrix3x3 RotationMatrix;
 protected:
 private:

    double s;                                                                   /// Real component
    dVector3 v;                                                                 /// Imaginary components
};
}
#endif