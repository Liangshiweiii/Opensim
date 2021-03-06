#include "Tools/Quaternion.h"

namespace opensim
{
    
void Quaternion::setRotationMatrix()
{
    if(length() < DBL_EPSILON)
    {
        RotationMatrix.set_to_unity();
    }    
    else
    {
        Quaternion Qcopy = normalized();

        RotationMatrix(0,0) = 1.0 - 2.0*(Qcopy.v[1]*Qcopy.v[1] + Qcopy.v[2]*Qcopy.v[2]);
        RotationMatrix(0,1) =       2.0*(Qcopy.v[0]*Qcopy.v[1] -    Qcopy.s*Qcopy.v[2]);
        RotationMatrix(0,2) =       2.0*(Qcopy.v[0]*Qcopy.v[2] +    Qcopy.s*Qcopy.v[1]);
        RotationMatrix(1,0) =       2.0*(Qcopy.v[0]*Qcopy.v[1] +    Qcopy.s*Qcopy.v[2]);
        RotationMatrix(1,1) = 1.0 - 2.0*(Qcopy.v[0]*Qcopy.v[0] + Qcopy.v[2]*Qcopy.v[2]);
        RotationMatrix(1,2) =       2.0*(Qcopy.v[1]*Qcopy.v[2] -    Qcopy.s*Qcopy.v[0]);
        RotationMatrix(2,0) =       2.0*(Qcopy.v[0]*Qcopy.v[2] -    Qcopy.s*Qcopy.v[1]);
        RotationMatrix(2,1) =       2.0*(Qcopy.v[1]*Qcopy.v[2] +    Qcopy.s*Qcopy.v[0]);
        RotationMatrix(2,2) = 1.0 - 2.0*(Qcopy.v[0]*Qcopy.v[0] + Qcopy.v[1]*Qcopy.v[1]);
    }
}

Quaternion Quaternion::lerp(const Quaternion& rhSQ1, const Quaternion& rhSQ2, const double t)
{
    // linear interpolation between two quaternions, where 0.0 < t < 1.0
    Quaternion result = (rhSQ1*(1.0 - t) + rhSQ2*t).normalized();
    return result;
}

Quaternion Quaternion::slerp(const Quaternion& Qa, const Quaternion& Qb, const double t)
{
    // linear interpolation between two quaternions, where 0.0 < t < 1.0
    // taken from http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/slerp/index.htm
    Quaternion result;

    double cosHalfTheta = Qa.s    * Qb.s 
                        + Qa.v[0] * Qb.v[0]
                        + Qa.v[1] * Qb.v[1]
                        + Qa.v[2] * Qb.v[2];
    double sign = 1.0;
    if (cosHalfTheta < 0)
    {
        sign = -1.0;
    }
    if(fabs(cosHalfTheta) >= 1.0)                                               // both rotations equal
    {
        result = Qa;
        return result;
    }
    else
    {
        double halfTheta = acos(sign*cosHalfTheta);
        double sinHalfTheta = sqrt(1.0 - cosHalfTheta*cosHalfTheta);
        if(fabs(sinHalfTheta) < 2.0*DBL_EPSILON)                                // the rotations are opposite of each other
        {
            result.s    = (Qa.s   *0.5 + sign*Qb.s   *0.5);
            result.v[0] = (Qa.v[0]*0.5 + sign*Qb.v[0]*0.5);
            result.v[1] = (Qa.v[1]*0.5 + sign*Qb.v[1]*0.5);
            result.v[2] = (Qa.v[2]*0.5 + sign*Qb.v[2]*0.5);
        }
        else                                                                    // general case
        {
            double ratioA = sin((1 - t) * halfTheta) / sinHalfTheta;
            double ratioB = sin(t * halfTheta) / sinHalfTheta;

            result.s    = (Qa.s    * ratioA + sign*Qb.s    * ratioB);
            result.v[0] = (Qa.v[0] * ratioA + sign*Qb.v[0] * ratioB);
            result.v[1] = (Qa.v[1] * ratioA + sign*Qb.v[1] * ratioB);
            result.v[2] = (Qa.v[2] * ratioA + sign*Qb.v[2] * ratioB);
        }
        result.setRotationMatrix();
    }
    return result;
}

std::string Quaternion::print(void) const
{
    std::stringstream out;
    out << "{ "   << s
        << ", < " << v[0]
        << ", "   << v[1]
        << ", "   << v[2] << " > }";
    return out.str();
}

std::string Quaternion::print_matrix(void) const
{
    std::stringstream out;
    out << RotationMatrix.print();
    return out.str();
}

}//namespace opensim