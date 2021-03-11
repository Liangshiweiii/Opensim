
#ifndef D3Q27_H
#define D3Q27_H

#include "Tools/Includes.h"

namespace opensim
{

class D3Q27                                                                     /// Lattice Boltzmann populations storage/manipulator
{
 public:

    static const double weights[3][3][3];

    double& operator()(int x, int y, int z)
    {
#ifdef DEBUG
        if(abs(x) > 1 or abs(y) > 1 or abs(z) > 1)
        {
            std::cout << "Error in D3Q27::operator()\n"
                 << "Access beyond storage range. (|x|, |y|, |z|) = ("
                 << abs(x) << "," << abs(y) << "," << abs(z) << ")"
                 << " > (1, 1, 1) "
                 << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[x + 1][y + 1][z + 1];
    };

    double operator()(int x, int y, int z) const
    {
#ifdef DEBUG
        if(abs(x) > 1 or abs(y) > 1 or abs(z) > 1)
        {
            std::cout << "Error in D3Q27::operator()\n"
                 << "Access beyond storage range. (|x|, |y|, |z|) = ("
                 << abs(x) << "," << abs(y) << "," << abs(z) << ")"
                 << " > (1, 1, 1) "
                 << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[x + 1][y + 1][z + 1];
    };

    D3Q27& operator=(D3Q27 rhs)
    {
        memmove(reinterpret_cast<double*> (storage), rhs.data(), 27*sizeof(double));
        return *this;
    };

    void set_to_zero()
    {
        memset(storage, 0, 27*sizeof(double));
    };

    D3Q27 operator* (const double value) const
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[1+x][1+y][1+z] * value;
        }
        return locPopulations;
    }

    D3Q27& operator*=(const double value)
    {
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            storage[1+x][1+y][1+z] *= value;
        }
        return *this;
    }

    D3Q27 operator/ (const double value) const
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[1+x][1+y][1+z] / value;
        }
        return locPopulations;
    }

    D3Q27& operator/=(const double value)
    {
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            storage[1+x][1+y][1+z] /= value;
        }
        return *this;
    }

    D3Q27 operator+ (const D3Q27 rhs) const
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[1+x][1+y][1+z] + rhs(x,y,z);
        }
        return locPopulations;
    };

    D3Q27& operator+=(const D3Q27 rhs)
    {
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            storage[1+x][1+y][1+z] += rhs(x,y,z);
        }
        return *this;
    };

    D3Q27 operator- (const D3Q27 rhs) const
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[1+x][1+y][1+z] - rhs(x,y,z);
        }
        return locPopulations;
    };

    D3Q27& operator-=(const D3Q27 rhs)
    {
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            storage[1+x][1+y][1+z] -= rhs(x,y,z);
        }
        return *this;
    };

    D3Q27 inverted(void)
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[1-x][1-y][1-z];
        }
        return locPopulations;
    };

    D3Q27 inverted(int x, int y, int z)
    {
        D3Q27 locPopulations = *this;
        locPopulations(x,y,z) = storage[1-x][1-y][1-z];

        return locPopulations;
    };

    D3Q27 Xreflected(void)
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[1-x][1+y][1+z];
        }
        return locPopulations;
    };

    D3Q27 Yreflected(void)
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[1+x][1-y][1+z];
        }
        return locPopulations;
    };

    D3Q27 Zreflected(void)
    {
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            locPopulations(x,y,z) = storage[1+x][1+y][1-z];
        }
        return locPopulations;
    };

    static D3Q27 EqDistribution(double lbDensity,
            dVector3 lbMomentum = {0.0,0.0,0.0})
    {
        dVector3 macroVel;
        if(lbDensity != 0)
        {
            macroVel = lbMomentum/lbDensity;
        }
        else
        {
            macroVel.set_to_zero();
        }
        double u2 = (macroVel*macroVel);
        D3Q27 locPopulations;
        for(int x = -1; x <= 1; x++)
        for(int y = -1; y <= 1; y++)
        for(int z = -1; z <= 1; z++)
        {
            double cu = x*macroVel[0] + y*macroVel[1] + z*macroVel[2];
            locPopulations(x,y,z) = lbDensity*weights[x+1][y+1][z+1]*(1.0 - 1.5*u2 + cu*(3.0 + 4.5*cu));
        }
        return locPopulations;
    };
    double* data()
    {
        return reinterpret_cast<double*>(storage);
    };

 protected:
 private:

    double storage[3][3][3];
};

} //namespace openphase
#endif
