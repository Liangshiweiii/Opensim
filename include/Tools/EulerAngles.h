#ifndef EULERANGLES_H
#define EULERANGLES_H

#include "Tools/Includes.h"

namespace opensim
{

enum EulerConvention{ XZX, XYX, YXY, YZY, ZXZ, ZYZ, /* proper Euler angles*/
                      XYZ, YZX, ZXY, XZY, ZYX, YXZ, /* Tait–Bryan angles*/
                      NNN /*default convention -> convention not set*/};

const std::vector<std::string>
     EulerConventionS{"XZX", "XYX", "YXY", "YZY", "ZXZ", "ZYZ", /* proper Euler angles*/
                      "XYZ", "YZX", "ZXY", "XZY", "ZYX", "YXZ", /* Tait–Bryan angles*/
                      "NNN" /*default convention -> convention not set*/};

class Quaternion;

class EulerAngles                                                               ///< Orientation angles and their Cos() and Sin().
{
//  Storages:
 public:
    union
    {
        double Q[3];
    };

    union
    {
        double CosQ[3];
    };

    union
    {
        double SinQ[3];
    };

    EulerConvention Convention;                                                 ///< Stores the convention for the Euler angles

    bool IsSet;                                                                 ///< Indicates whether sin/cos are set.

//  Methods:

    EulerAngles(): IsSet(false)                                                 ///< Constructor
    {
        Convention = NNN;
        
        Q[0] = 0.0;
        Q[1] = 0.0;
        Q[2] = 0.0;

        SinQ[0] = 0.0;
        SinQ[1] = 0.0;
        SinQ[2] = 0.0;

        CosQ[0] = 1.0;
        CosQ[1] = 1.0;
        CosQ[2] = 1.0;
    };

    EulerAngles(std::initializer_list<double> Angles, EulerConvention locConvention)
    {
#ifdef DEBUG
        if (Angles.size() != 3)
        {
            std::cout << "Error in EulerAngles::constructor()\n"
                      << "Initialization list size beyond storage range."
                      << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif

        unsigned int ii = 0;
        for (auto it = Angles.begin(); it != Angles.end(); it++)
        {
            Q[ii] = *it;
            ii += 1;
        }
        
        Convention = locConvention;

        SinQ[0] = sin(Q[0]);
        SinQ[1] = sin(Q[1]);
        SinQ[2] = sin(Q[2]);

        CosQ[0] = cos(Q[0]);
        CosQ[1] = cos(Q[1]);
        CosQ[2] = cos(Q[2]);

        IsSet = true;
    }

    EulerAngles(const EulerAngles& rhs)
    {
        Convention = rhs.Convention;

        Q[0] = rhs.Q[0];
        Q[1] = rhs.Q[1];
        Q[2] = rhs.Q[2];

        SinQ[0] = rhs.SinQ[0];
        SinQ[1] = rhs.SinQ[1];
        SinQ[2] = rhs.SinQ[2];

        CosQ[0] = rhs.CosQ[0];
        CosQ[1] = rhs.CosQ[1];
        CosQ[2] = rhs.CosQ[2];

        IsSet = rhs.IsSet;
    };

    bool operator==(const EulerAngles& rhs)
    {
        double epsilon = fabs(Q[0] - rhs.Q[0]) + fabs(Q[1] - rhs.Q[1]) + fabs(Q[2] - rhs.Q[2]);
        if(Convention == rhs.Convention and epsilon < 3.0*DBL_EPSILON)
        {
            return true;
        }
        return false;
    }

    void set(const double q1, const double q2, const double q3, EulerConvention locConvention)
    {
        Convention = locConvention;

        Q[0] = q1;
        Q[1] = q2;
        Q[2] = q3;

        SinQ[0] = sin(Q[0]);
        SinQ[1] = sin(Q[1]);
        SinQ[2] = sin(Q[2]);

        CosQ[0] = cos(Q[0]);
        CosQ[1] = cos(Q[1]);
        CosQ[2] = cos(Q[2]);

        IsSet = true;
    };
    
    void setTrigonometricFunctions()
    {
        if(IsSet == false)
        {
            SinQ[0] = sin(Q[0]);
            SinQ[1] = sin(Q[1]);
            SinQ[2] = sin(Q[2]);

            CosQ[0] = cos(Q[0]);
            CosQ[1] = cos(Q[1]);
            CosQ[2] = cos(Q[2]);

            IsSet = true;
        }
    };

    void set_to_zero(void)
    {
        if(Convention != NNN)
        {
            Q[0] = 0.0;
            Q[1] = 0.0;
            Q[2] = 0.0;

            SinQ[0] = 0.0;
            SinQ[1] = 0.0;
            SinQ[2] = 0.0;

            CosQ[0] = 1.0;
            CosQ[1] = 1.0;
            CosQ[2] = 1.0;

            IsSet = false;
        }
        else
        {
            std::cout << " EulerAngles::set_to_zero(): Trying to set values of"
                      << " EulerAngles object that has no valid convention!"
                      << " Use EulerAngles::set_to_zero(q1, q2, q3, Convention)"
                      << " instead! Terminating!" << std::endl;
            exit(13);
        }
    };

    void set_convention(const EulerConvention locConvention)
    {
        if(Convention == NNN)
        {
            Convention = locConvention;
        }
        else
        {
            std::cout << " EulerAngles::set_convention(): Trying to set convention of"
                      << " EulerAngles object that has a valid convention already!"
                      << " Use EulerAngles::set(q1, q2, q3, Convention)"
                      << " instead! Terminating!"<< std::endl;
            exit(13);

        }
    };

    void add(const double q1, const double q2, const double q3)
    {
        Q[0] += q1;
        Q[1] += q2;
        Q[2] += q3;

        SinQ[0] = sin(Q[0]);
        SinQ[1] = sin(Q[1]);
        SinQ[2] = sin(Q[2]);

        CosQ[0] = cos(Q[0]);
        CosQ[1] = cos(Q[1]);
        CosQ[2] = cos(Q[2]);

        IsSet = true;
    };

    EulerAngles& operator=(const EulerAngles& rhs)
    {
#ifdef DEBUG
        if (Convention != rhs.Convention and 
            Convention != NNN)
        {
            std::cout << "Euler angle conventions do not coincide.\n"
                      << "EulerAngles::operator=.\n"
                      << "Terminating!" << std::endl;
            exit(1);
        }
#endif

        Convention = rhs.Convention;

        Q[0] = rhs.Q[0];
        Q[1] = rhs.Q[1];
        Q[2] = rhs.Q[2];

        SinQ[0] = sin(rhs.Q[0]);
        SinQ[1] = sin(rhs.Q[1]);
        SinQ[2] = sin(rhs.Q[2]);

        CosQ[0] = cos(rhs.Q[0]);
        CosQ[1] = cos(rhs.Q[1]);
        CosQ[2] = cos(rhs.Q[2]);

        IsSet = true;

        return *this;
    };

    EulerAngles operator+(const EulerAngles& rhs) const
    {
#ifdef DEBUG
        if (Convention != rhs.Convention)
        {
            std::cout << "Euler angle convention do not coincide. Method EulerAngles::operator+. Exit!" << std::endl;
            exit(1);
        }
#endif
        EulerAngles returnAng;

        returnAng.Convention = rhs.Convention;
        
        returnAng.Q[0] = Q[0] + rhs.Q[0];
        returnAng.Q[1] = Q[1] + rhs.Q[1];
        returnAng.Q[2] = Q[2] + rhs.Q[2];

        returnAng.SinQ[0] = sin(returnAng.Q[0]);
        returnAng.SinQ[1] = sin(returnAng.Q[1]);
        returnAng.SinQ[2] = sin(returnAng.Q[2]);

        returnAng.CosQ[0] = cos(returnAng.Q[0]);
        returnAng.CosQ[1] = cos(returnAng.Q[1]);
        returnAng.CosQ[2] = cos(returnAng.Q[2]);

        returnAng.IsSet = true;

        return returnAng;
    };

    EulerAngles operator-(const EulerAngles& rhs) const
    {
#ifdef DEBUG
        if (Convention != rhs.Convention)
        {
            std::cout << "Euler angle convention do not coincide. Method EulerAngles::operator-. Exit!" << std::endl;
            exit(1);
        }
#endif
        EulerAngles returnAng;

        returnAng.Convention = rhs.Convention;
        
        returnAng.Q[0] = Q[0] - rhs.Q[0];
        returnAng.Q[1] = Q[1] - rhs.Q[1];
        returnAng.Q[2] = Q[2] - rhs.Q[2];

        returnAng.SinQ[0] = sin(returnAng.Q[0]);
        returnAng.SinQ[1] = sin(returnAng.Q[1]);
        returnAng.SinQ[2] = sin(returnAng.Q[2]);

        returnAng.CosQ[0] = cos(returnAng.Q[0]);
        returnAng.CosQ[1] = cos(returnAng.Q[1]);
        returnAng.CosQ[2] = cos(returnAng.Q[2]);

        returnAng.IsSet = true;

        return returnAng;
    };

    EulerAngles& operator+=(const EulerAngles& rhs)
    {
#ifdef DEBUG
        if (Convention != rhs.Convention)
        {
            std::cout << "Euler angle convention do not coincide. Method EulerAngles::operator +=. Exit!" << std::endl;
            exit(1);
        }
#endif
        Q[0] += rhs.Q[0];
        Q[1] += rhs.Q[1];
        Q[2] += rhs.Q[2];

        SinQ[0] = sin(Q[0]);
        SinQ[1] = sin(Q[1]);
        SinQ[2] = sin(Q[2]);

        CosQ[0] = cos(Q[0]);
        CosQ[1] = cos(Q[1]);
        CosQ[2] = cos(Q[2]);

        IsSet = true;

        return *this;
    };

    EulerAngles& operator-=(const EulerAngles& rhs)
    {
#ifdef DEBUG
        if (Convention != rhs.Convention)
        {
            std::cout << "Euler angle convention do not coincide. Method EulerAngles::operator-=. Exit!" << std::endl;
            exit(1);
        }
#endif
        Q[0] += rhs.Q[0];
        Q[1] += rhs.Q[1];
        Q[2] += rhs.Q[2];

        SinQ[0] = sin(Q[0]);
        SinQ[1] = sin(Q[1]);
        SinQ[2] = sin(Q[2]);

        CosQ[0] = cos(Q[0]);
        CosQ[1] = cos(Q[1]);
        CosQ[2] = cos(Q[2]);

        IsSet = true;

        return *this;
    };

    EulerAngles get_degree(void) const
    {
        EulerAngles returnAng;

        returnAng.Convention = Convention;

        returnAng.SinQ[0] = sin(Q[0]);
        returnAng.SinQ[1] = sin(Q[1]);
        returnAng.SinQ[2] = sin(Q[2]);

        returnAng.CosQ[0] = cos(Q[0]);
        returnAng.CosQ[1] = cos(Q[1]);
        returnAng.CosQ[2] = cos(Q[2]);

        returnAng.Q[0] = Q[0]*180.0/3.14159265358979323846;
        returnAng.Q[1] = Q[1]*180.0/3.14159265358979323846;
        returnAng.Q[2] = Q[2]*180.0/3.14159265358979323846;

        returnAng.IsSet = true;

        return returnAng;
    };

    std::string get_convention(void) const
    {
        std::string returnCon = EulerConventionS[Convention];
        return returnCon;
    };

    EulerAngles operator*(const double rhs) const
    {
        EulerAngles returnAng;
        returnAng.set_convention(Convention);
        returnAng.set_to_zero();

        returnAng.Q[0] = Q[0]*rhs;
        returnAng.Q[1] = Q[1]*rhs;
        returnAng.Q[2] = Q[2]*rhs;

        returnAng.SinQ[0] = sin(returnAng.Q[0]);
        returnAng.SinQ[1] = sin(returnAng.Q[1]);
        returnAng.SinQ[2] = sin(returnAng.Q[2]);

        returnAng.CosQ[0] = cos(returnAng.Q[0]);
        returnAng.CosQ[1] = cos(returnAng.Q[1]);
        returnAng.CosQ[2] = cos(returnAng.Q[2]);

        return returnAng;
    };

    EulerAngles& operator*=(const double rhs)
    {
        Q[0] = Q[0]*rhs;
        Q[1] = Q[1]*rhs;
        Q[2] = Q[2]*rhs;

        SinQ[0] = sin(Q[0]);
        SinQ[1] = sin(Q[1]);
        SinQ[2] = sin(Q[2]);

        CosQ[0] = cos(Q[0]);
        CosQ[1] = cos(Q[1]);
        CosQ[2] = cos(Q[2]);

        return *this;
    };
    
    dMatrix3x3 getRotationMatrix() const;
    Quaternion getQuaternion() const;
    
    std::string print(void) const
    {
        std::stringstream out;

        out << "[" << Q[0] << ", "
                   << Q[1] << ", "
                   << Q[2] << "]["
                   << EulerConventionS[Convention]
                   << "]";
       return out.str();
    };
    std::string print_degree(void) const
    {
        std::stringstream out;

        out << "[" << Q[0]* 180.0/3.14159265358979323846 << ", "
                   << Q[1]* 180.0/3.14159265358979323846 << ", "
                   << Q[2]* 180.0/3.14159265358979323846 << "]["
                   << EulerConventionS[Convention]
                   << "]";
       return out.str();
    };
    std::string print_entire(void) const
    {
        std::stringstream out;

        out << "Angle      [" << Q[0] << ", "
                              << Q[1] << ", "
                              << Q[2] << "]" << std::endl
            << "Convention [" << EulerConventionS[Convention]
                              << "]" << std::endl
            << "Sin        [" << SinQ[0] << ", "
                              << SinQ[1] << ", "
                              << SinQ[2] << "]" << std::endl
            << "Cos        [" << CosQ[0] << ", "
                              << CosQ[1] << ", "
                              << CosQ[2] << "]" << std::endl
            << "Is set:     " << IsSet << std::endl;
       return out.str();
    };
 protected:
 private:
};

}// namespace opensim
#endif
