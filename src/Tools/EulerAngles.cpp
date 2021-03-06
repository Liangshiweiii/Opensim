#include "Tools/Includes.h"
#include "Tools/EulerAngles.h"
#include "Tools/Quaternion.h"
#include "Info.h"

namespace opensim
{
    
dMatrix3x3 EulerAngles::getRotationMatrix() const
{
    dMatrix3x3 result;
    double c1 = CosQ[0];
    double c2 = CosQ[1];
    double c3 = CosQ[2];
    double s1 = SinQ[0];
    double s2 = SinQ[1];
    double s3 = SinQ[2];

    //The rotation matrices follow notations from http://en.wikipedia.org/wiki/Euler_angles
    switch (Convention)
    {
        // Proper Euler angles
        case XZX:
        {
            double XZX[3][3] = 
            {{                 c2,            - c3*s2,              s2*s3},
             {              c1*s2,   c1*c2*c3 - s1*s3, - c3*s1 - c1*c2*s3},
             {              s1*s2,   c1*s3 + c2*c3*s1,   c1*c3 - c2*s1*s3}};
             
            memmove(result.data(), XZX, 9*sizeof(double));
            break;
        }
        case XYX:
        {
            double XYX[3][3] = 
            {{                 c2,             s2*s3,              c3*s2},
             {              s1*s2,   c1*c3 - c2*s1*s3, - c1*s3 - c2*c3*s1},
             {            - c1*s2,   c3*s1 + c1*c2*s3,   c1*c2*c3 - s1*s3}};
             
            memmove(result.data(), XYX, 9*sizeof(double));
            break;
        }
        case YXY:
        {
            double YXY[3][3] = 
            {{   c1*c3 - c2*s1*s3,              s1*s2,   c2*c3*s1 + c1*s3},
             {              s2*s3,                 c2,            - c3*s2},
             { - c3*s1 - c1*c2*s3,              c1*s2,   c1*c2*c3 - s1*s3}};
            memmove(result.data(), YXY, 9*sizeof(double));
            break;
        }
        case YZY:
        {
            double YZY[3][3] = 
            {{   c1*c2*c3 - s1*s3,            - c1*s2,   c3*s1 + c1*c2*s3},
             {              c3*s2,                 c2,              s2*s3},
             {-(c2*c3*s1) - c1*s3,              s1*s2,   c1*c3 - c2*s1*s3}};
            memmove(result.data(), YZY, 9*sizeof(double));
            break;
        }
        case ZYZ:
        {
            double ZYZ[3][3] =
            {{   c1*c2*c3 - s1*s3, - c3*s1 - c1*c2*s3,             c1*s2},
             {   c2*c3*s1 + c1*s3,   c1*c3 - c2*s1*s3,             s1*s2},
             {            - c3*s2,              s2*s3,                c2}};
            memmove(result.data(), ZYZ, 9*sizeof(double));
            break;
        }
        case ZXZ:
        {
            double ZXZ[3][3] =
            {{   c1*c3 - c2*s1*s3, - c2*c3*s1 - c1*s3,             s1*s2},
             {   c3*s1 + c1*c2*s3,   c1*c2*c3 - s1*s3,           - c1*s2},
             {              s2*s3,              c3*s2,                c2}};
            memmove(result.data(), ZXZ, 9*sizeof(double));
            break;
        }     
                
        // Tait-Bryan angles
        case XYZ:
        {
            double XYZ[3][3] =
            {{              c2*c3,            - c2*s3,                s2},
             {   c1*s3 + c3*s1*s2,   c1*c3 - s1*s2*s3,           - c2*s1},
             {   s1*s3 - c1*c3*s2,   c1*s2*s3 + c3*s1,             c1*c2}};
            memmove(result.data(), XYZ, 9*sizeof(double));
            break;
        }
        case ZYX:
        {
            double ZYX[3][3] =
            {{              c1*c2,  c1*s2*s3 - c3*s1,   s1*s3 + c1*c3*s2},
             {              c2*s1,  c1*c3 - s1*s2*s3,   c3*s1*s2 - c1*s3},
             {               - s2,             c2*s3,              c2*c3}};
            memmove(result.data(), ZYX, 9*sizeof(double));
            break;
        }
        default:
        {
            std::cout << " EulerAngles::getRotationMatrix(): Wrong Euler convention is used!"
              << " Terminating!"<< std::endl;
            exit(13);        
        }
    }    
    return result;
}

Quaternion EulerAngles::getQuaternion(void) const
{
    // Conversions taken from the following NASA document:
    // http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19770024290.pdf
    // Flipped first and last angle to have active rotation quaternion instead
    // of passive rotation as described in reference above.

    double Eang1  = Q[2];
    //double Eang2 = Q[1];
    double Eang3  = Q[0];
    double Eang1h = 0.5*Q[2];
    double Eang2h = 0.5*Q[1];
    double Eang3h = 0.5*Q[0];

    Quaternion result;

    switch (Convention)
    {
        case ZYX:
        {
            result.set(-sin(Eang1h)*sin(Eang2h)*sin(Eang3h) + 
                        cos(Eang1h)*cos(Eang2h)*cos(Eang3h),
                        sin(Eang1h)*cos(Eang2h)*cos(Eang3h) + 
                        sin(Eang2h)*sin(Eang3h)*cos(Eang1h),
                       -sin(Eang1h)*sin(Eang3h)*cos(Eang2h) + 
                        sin(Eang2h)*cos(Eang1h)*cos(Eang3h),
                        sin(Eang1h)*sin(Eang2h)*cos(Eang3h) + 
                        sin(Eang3h)*cos(Eang1h)*cos(Eang2h));
            break;
        }
        case YZX:
        {
            result.set( sin(Eang1h)*sin(Eang2h)*sin(Eang3h) + 
                        cos(Eang1h)*cos(Eang2h)*cos(Eang3h),
                        sin(Eang1h)*cos(Eang2h)*cos(Eang3h) - 
                        sin(Eang2h)*sin(Eang3h)*cos(Eang1h),
                       -sin(Eang1h)*sin(Eang2h)*cos(Eang3h) + 
                        sin(Eang3h)*cos(Eang1h)*cos(Eang2h),
                        sin(Eang1h)*sin(Eang3h)*cos(Eang2h) + 
                        sin(Eang2h)*cos(Eang1h)*cos(Eang3h));
            break;
        }
        case XYX:
        {
            result.set(cos(Eang2h)*cos((Eang1 + Eang3)/2.0),
                       cos(Eang2h)*sin((Eang1 + Eang3)/2.0),
                       sin(Eang2h)*cos((Eang1 - Eang3)/2.0),
                       sin(Eang2h)*sin((Eang1 - Eang3)/2.0));
            break;
        }
        case XZX:
        {
            result.set(cos(Eang2h)*cos((Eang1 + Eang3)/2.0),
                       cos(Eang2h)*sin((Eang1 + Eang3)/2.0),
                      -sin(Eang2h)*sin((Eang1 - Eang3)/2.0),
                       sin(Eang2h)*cos((Eang1 - Eang3)/2.0));
            break;
        }
        case ZXY:
        {
            result.set( sin(Eang1h)*sin(Eang2h)*sin(Eang3h) + 
                        cos(Eang1h)*cos(Eang2h)*cos(Eang3h),
                        sin(Eang1h)*sin(Eang3h)*cos(Eang2h) + 
                        sin(Eang2h)*cos(Eang1h)*cos(Eang3h),
                        sin(Eang1h)*cos(Eang2h)*cos(Eang3h) - 
                        sin(Eang2h)*sin(Eang3h)*cos(Eang1h),
                       -sin(Eang1h)*sin(Eang2h)*cos(Eang3h) + 
                        sin(Eang3h)*cos(Eang1h)*cos(Eang2h));
            break;
        }
        case XZY:
        {
            result.set(-sin(Eang1h)*sin(Eang2h)*sin(Eang3h) + 
                        cos(Eang1h)*cos(Eang2h)*cos(Eang3h),
                        sin(Eang1h)*sin(Eang2h)*cos(Eang3h) + 
                        sin(Eang3h)*cos(Eang1h)*cos(Eang2h),
                        sin(Eang1h)*cos(Eang2h)*cos(Eang3h) + 
                        sin(Eang2h)*sin(Eang3h)*cos(Eang1h),
                       -sin(Eang1h)*sin(Eang3h)*cos(Eang2h) + 
                        sin(Eang2h)*cos(Eang1h)*cos(Eang3h));
            break;
        }
        case YXY:
        {
            result.set(cos(Eang2h)*cos((Eang1+Eang3)/2.0),
                       sin(Eang2h)*cos((Eang1-Eang3)/2.0),
                       cos(Eang2h)*sin((Eang1+Eang3)/2.0),
                      -sin(Eang2h)*sin((Eang1-Eang3)/2.0));
            break;
        }
        case YZY:
        {
            result.set(cos(Eang2h)*cos((Eang1+Eang3)/2.0),
                       sin(Eang2h)*sin((Eang1-Eang3)/2.0),
                       cos(Eang2h)*sin((Eang1+Eang3)/2.0),
                       sin(Eang2h)*cos((Eang1-Eang3)/2.0));
            break;
        }
        case YXZ:
        {
            result.set(-sin(Eang1h)*sin(Eang2h)*sin(Eang3h) + 
                        cos(Eang1h)*cos(Eang2h)*cos(Eang3h),
                       -sin(Eang1h)*sin(Eang2h)*cos(Eang3h) + 
                        sin(Eang2h)*cos(Eang1h)*cos(Eang3h),
                        sin(Eang1h)*sin(Eang2h)*cos(Eang3h) + 
                        sin(Eang3h)*cos(Eang1h)*cos(Eang2h),
                        sin(Eang1h)*cos(Eang2h)*cos(Eang3h) + 
                        sin(Eang2h)*sin(Eang3h)*cos(Eang1h));
            break;
        }
        case XYZ:
        {
            result.set( sin(Eang1h)*sin(Eang2h)*sin(Eang3h) + 
                        cos(Eang1h)*cos(Eang2h)*cos(Eang3h),
                       -sin(Eang1h)*sin(Eang2h)*cos(Eang3h) + 
                        sin(Eang3h)*cos(Eang1h)*cos(Eang2h),
                        sin(Eang1h)*sin(Eang3h)*cos(Eang2h) + 
                        sin(Eang2h)*cos(Eang1h)*cos(Eang3h),
                        sin(Eang1h)*cos(Eang2h)*cos(Eang3h) + 
                        sin(Eang2h)*sin(Eang3h)*cos(Eang1h));
            break;
        }
        case ZXZ:
       {
            result.set(cos(Eang2h)*cos((Eang1+Eang3)/2.0),
                       sin(Eang2h)*cos((Eang1-Eang3)/2.0),
                       sin(Eang2h)*sin((Eang1-Eang3)/2.0),
                       cos(Eang2h)*sin((Eang1+Eang3)/2.0));
            break;
        }
        case ZYZ:
        {
            result.set(cos(Eang2h)*cos((Eang1+Eang3)/2.0),
                      -sin(Eang2h)*sin((Eang1-Eang3)/2.0),
                       sin(Eang2h)*cos((Eang1-Eang3)/2.0),
                       cos(Eang2h)*sin((Eang1+Eang3)/2.0));
            break;
        }
        default:
        {
            Info::WriteExit("Wrong or not supported convention.",
                            "EulerAngles", "getQuaternion()");
            exit(1);
        }
    }
    return result;
}

}//namespace openphase
