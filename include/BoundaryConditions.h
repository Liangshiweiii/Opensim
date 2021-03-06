#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H

#include "Tools/Includes.h"

namespace opensim
{
    class Settings;
    /// Types of boundary conditions.

    const long int Periodic = 0;
    const long int NoFlux   = 1;
    const long int Free     = 2;
    const long int Fixed    = 3;

    class BoundaryConditions : public OPObject
    {
        public:

        BoundaryConditions(){};
        BoundaryConditions(const Settings& locSettings, const std::string InputFileName = "NONE");

        template< class T, long int Num >
        void SetX(Storage3D<T, Num>& loc3Dstorage) const;                           /// 设置沿X边界存储标准值的边界条件
        template< class T, long int Num >
        void SetY(Storage3D<T, Num>& loc3Dstorage) const;                           /// 设置沿Y边界存储标准值的边界条件
        template< class T, long int Num >
        void SetZ(Storage3D<T, Num>& loc3Dstorage) const;                           /// 设置沿Z边界存储标准值的边界条件

        template< class T>
        void SetXVector(Storage3D<T, 0>& loc3Dstorage) const;                       /// 设置沿X边界存储标准值的边界条件
        template< class T>
        void SetYVector(Storage3D<T, 0>& loc3Dstorage) const;                       /// 设置沿Y边界存储标准值的边界条件
        template< class T>
        void SetZVector(Storage3D<T, 0>& loc3Dstorage) const;                       /// 设置沿Z边界存储标准值的边界条件

    /*    template< class T>
        void SetXVectorNoFlux(Storage3D<T, 0>& loc3Dstorage) const;                 /// 设置沿X边界存储标准值的边界条件
        template< class T>
        void SetYVectorNoFlux(Storage3D<T, 0>& loc3Dstorage) const;                 /// 设置沿Y边界存储标准值的边界条件
        template< class T>
        void SetZVectorNoFlux(Storage3D<T, 0>& loc3Dstorage) const;                 /// 设置沿Z边界存储标准值的边界条件

        template< class T>
        void SetLBXVector(Storage3D<T, 0>& loc3Dstorage) const;                     /// 设置沿X边界存储标准值的边界条件ies
        template< class T>
        void SetLBYVector(Storage3D<T, 0>& loc3Dstorage) const;                     /// 设置沿Y边界存储标准值的边界条件ies
        template< class T>
        void SetLBZVector(Storage3D<T, 0>& loc3Dstorage) const;                     /// 设置沿Z边界存储标准值的边界条件ies
    */
        template< class T>
        void SetXVector(Storage3D<T, 1>& loc3Dstorage) const;                       /// 设置沿X边界存储标准值的边界条件
        template< class T>
        void SetYVector(Storage3D<T, 1>& loc3Dstorage) const;                       /// 设置沿Y边界存储标准值的边界条件
        template< class T>
        void SetZVector(Storage3D<T, 1>& loc3Dstorage) const;                       /// 设置沿Z边界存储标准值的边界条件

        template< class T>
        void SetXVector(Storage3D<T, 2>& loc3Dstorage) const;                       /// 设置沿X边界存储标准值的边界条件
        template< class T>
        void SetYVector(Storage3D<T, 2>& loc3Dstorage) const;                       /// 设置沿Y边界存储标准值的边界条件
        template< class T>
        void SetZVector(Storage3D<T, 2>& loc3Dstorage) const;                       /// 设置沿Z边界存储标准值的边界条件

        template< class T>
        void SetXVector(Storage3D<T, 3>& loc3Dstorage) const;                       /// 设置沿X边界存储标准值的边界条件
        template< class T>
        void SetYVector(Storage3D<T, 3>& loc3Dstorage) const;                       /// 设置沿Y边界存储标准值的边界条件
        template< class T>
        void SetZVector(Storage3D<T, 3>& loc3Dstorage) const;                       /// 设置沿Z边界存储标准值的边界条件

        template< class T, long int Num >
        void SetXFlags(Storage3D<T, Num>& loc3Dstorage) const;                      /// 设置边界条件以沿X边界存储标志
        template< class T, long int Num >
        void SetYFlags(Storage3D<T, Num>& loc3Dstorage) const;                      /// 设置边界条件以沿Y边界存储标志
        template< class T, long int Num >
        void SetZFlags(Storage3D<T, Num>& loc3Dstorage) const;                      /// 设置边界条件以沿Z边界存储标志
        using OPObject::ReadInput;
        void Initialize(const Settings& Settings);                                  /// 初始化类的变量
        void ReadInput(std::string InputFileName);                                  /// 读取方向属性

        long int Index(const long int x, const long int y, const long int z,
                    const long int Nx, const long int Ny, const long int Nz,
                    const long int Bcells) const;                                /// 使用边界条件的（x，y，z）位置的索引。
        long int BC0X;
        long int BCNX;
        long int BC0Y;
        long int BCNY;
        long int BC0Z;
        long int BCNZ;

        protected:
        private:
    };

    inline long int BoundaryConditions::Index(const long int x, const long int y, const long int z,
                                        const long int Nx, const long int Ny, const long int Nz, const long int Bcells) const
        {
            long int xx = x;
            long int yy = y;
            long int zz = z;
            if(x < 0)
            {
                switch (BC0X)
                {
                    case NoFlux:
                    {
                        xx = (-x - 1);
                        break;

                    }
                    case Periodic:
                    {
                        xx = (x + Nx)%Nx;
                        break;
                    }
                    case Free:
                    case Fixed:
                    {
                        xx = 0;
                        break;
                    }
                    default:
                    {
                        xx = std::max(xx,-Bcells);
                        break;
                    }
                }
            }
            if(x > Nx - 1)
            {
                switch (BCNX)
                {
                    case NoFlux:
                    {
                        xx = 2*Nx - 1 - x;
                        break;
                    }
                    case Periodic:
                    {
                        xx = x%Nx;
                        break;
                    }
                    case Free:
                    case Fixed:
                    {
                        xx = Nx - 1;
                        break;
                    }
                    default:
                    {
                        xx = std::min(xx,Nx+Bcells-1);
                        break;
                    }
                }
            }

            if(y < 0)
            {
                switch (BC0Y)
                {
                    case NoFlux:
                    {
                        yy = (-y - 1);
                        break;

                    }
                    case Periodic:
                    {
                        yy = (y + Ny)%Ny;
                        break;
                    }
                    case Free:
                    case Fixed:
                    {
                        yy = 0;
                        break;
                    }
                    default:
                    {
                        yy = std::max(yy,-Bcells);
                        break;
                    }
                }
            }
            if(y > Ny - 1)
            {
                switch (BCNY)
                {
                    case NoFlux:
                    {
                        yy = 2*Ny - 1 - y;
                        break;
                    }
                    case Periodic:
                    {
                        yy = y%Ny;
                        break;
                    }
                    case Free:
                    case Fixed:
                    {
                        yy = Ny - 1;
                        break;
                    }
                    default:
                    {
                        yy = std::min(yy,Ny+Bcells-1);
                        break;
                    }
                }
            }

            if(z < 0)
            {
                switch (BC0Z)
                {
                    case NoFlux:
                    {
                        zz = (-z - 1);
                        break;

                    }
                    case Periodic:
                    {
                        zz = (z + Nz)%Nz;
                        break;
                    }
                    case Free:
                    case Fixed:
                    {
                        zz = 0;
                        break;
                    }
                    default:
                    {
                        zz = std::max(zz,-Bcells);
                        break;
                    }
                }
            }
            if(z > Nz - 1)
            {
                switch (BCNZ)
                {
                    case NoFlux:
                    {
                        zz = 2*Nz - 1 - z;
                        break;
                    }
                    case Periodic:
                    {
                        zz = z%Nz;
                        break;
                    }
                    case Free:
                    case Fixed:
                    {
                        zz = Nz - 1;
                        break;
                    }
                    default:
                    {
                        zz = std::min(zz,Nz+Bcells-1);
                        break;
                    }
                }
            }
            return ((Ny + 2*Bcells)*(xx + Bcells) + yy + Bcells)*(Nz + 2*Bcells) + zz + Bcells;
        }

    // Flags
    template< class T, long int Num >
    void BoundaryConditions::SetXFlags(Storage3D<T, Num> &Field) const
    {
        if(BC0X == Periodic || BCNX == Periodic)
        {
            for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
            for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
            {
                Field( 0, j, k) = std::max(Field( 0, j, k), Field(Field.sizeX(), j, k));
                Field(Field.sizeX() - 1, j, k) = std::max(Field(Field.sizeX() - 1, j, k), Field(-1, j, k));
            }

            for(long int i = -Field.Bcells(); i < 0; i++)
            for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
            for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
            {
                Field( i, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                Field(Field.sizeX() - i - 1, j, k) = Field((- i - 1)%Field.sizeX(), j, k);
            }
            return;
        }

        switch (BC0X)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field((- i - 1)%Field.sizeZ(), j, k);
                }
                break;
            }
            case Free:
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNX)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(Field.sizeX() - 1 - i, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                }
                break;
            }
            case Free:
            case Fixed:
            default:
            {
                break;
            }
        }
    }
    template< class T, long int Num >
    void BoundaryConditions::SetYFlags(Storage3D<T, Num> &Field) const
    {
        if(BC0Y == Periodic || BCNY == Periodic)
        {
            for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
            for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
            {
                Field(i, 0, k) = std::max(Field(i, 0, k), Field(i, Field.sizeY(), k));
                Field(i, Field.sizeY() - 1, k) = std::max(Field(i, Field.sizeY() - 1, k), Field(i, -1, k));
            }

            for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
            for(long int j = -Field.Bcells(); j < 0; j++)
            for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
            {
                Field(i, j, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                Field(i, Field.sizeY() - j -1, k) = Field(i, (- j - 1)%Field.sizeY(), k);
            }
            return;
        }

        switch (BC0Y)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field(i, (- j - 1)%Field.sizeY(), k);
                }
                break;
            }
            case Free:
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNY)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, Field.sizeY() - 1 - j, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                }
                break;
            }
            case Free:
            case Fixed:
            default:
            {
                break;
            }
        }
    }
    template< class T, long int Num >
    void BoundaryConditions::SetZFlags(Storage3D<T, Num> &Field) const
    {
        if(BC0Z == Periodic || BCNZ == Periodic)
        {
            for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
            for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
            {
                Field(i, j, 0) = std::max(Field(i, j, 0), Field(i, j, Field.sizeZ()));
                Field(i, j, Field.sizeZ() - 1) = std::max(Field(i, j, Field.sizeZ() - 1), Field(i, j, -1));
            }

            for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
            for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
            for(long int k = -Field.Bcells(); k < 0; k++)
            {
                Field(i, j, k) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (- k - 1)%Field.sizeZ());
            }
            return;
        }

        switch (BC0Z)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, (- k - 1)%Field.sizeZ());
                }
                break;
            }
            case Free:
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNZ)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, Field.sizeZ() - 1 - k) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                }
                break;
            }
            case Free:
            case Fixed:
            default:
            {
                break;
            }
        }
    }

    // Standard
    template< class T, long int Num >
    void BoundaryConditions::SetX(Storage3D<T, Num> &Field) const
    {
        if(BC0X == Periodic || BCNX == Periodic)
        {
            for(long int i = -Field.Bcells(); i < 0; i++)
            for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
            for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
            {
                Field(i, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                Field(Field.sizeX() - i - 1, j, k) = Field((- i - 1)%Field.sizeX(), j, k);
            }
            return;
        }

        switch (BC0X)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field((-i - 1)%Field.sizeX(), j, k);
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field(0, j, k)
                        + (Field(1%Field.sizeX(), j, k) - Field(0, j, k)) * i;
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNX)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(Field.sizeX() - i - 1, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    const long int i0 = Field.sizeX() - 1;
                    Field(i0 - i, j, k) = Field(i0, j, k)
                        +(Field(i0, j, k) - Field((i0-1)%Field.sizeX(), j, k))*(-i);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }
    template< class T, long int Num >
    void BoundaryConditions::SetY(Storage3D<T, Num> &Field) const
    {
        if(BC0Y == Periodic || BCNY == Periodic)
        {
            for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
            for(long int j = -Field.Bcells(); j < 0; j++)
            for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
            {
                Field(i, j, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                Field(i, Field.sizeY() - j - 1, k) = Field(i, (-j - 1)%Field.sizeY(), k);
            }
            return;
        }

        switch (BC0Y)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field(i, (-j - 1)%Field.sizeY(), k);
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field(i, 0, k)
                        + (Field(i, 1%Field.sizeY(), k) - Field(i, 0, k)) * j;
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNY)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, Field.sizeY() - j - 1, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    const long int j0 = Field.sizeY() - 1;
                    Field(i, j0 - j, k) = Field(i, j0, k)
                        +(Field(i, j0, k) - Field(i, (j0-1)%Field.sizeY(), k))*(-j);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }
    template< class T, long int Num >
    void BoundaryConditions::SetZ(Storage3D<T, Num> &Field) const
    {
        if(BC0Z == Periodic || BCNZ == Periodic)
        {
            for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
            for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
            for(long int k = -Field.Bcells(); k < 0; k++)
            {
                Field(i, j, k) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (-k - 1)%Field.sizeZ());
            }
            return;
        }

        switch (BC0Z)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, (-k - 1)%Field.sizeZ());
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, 0)
                        + (Field(i, j, 1%Field.sizeZ()) - Field(i, j, 0)) * k;
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNZ)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    const long int k0 = Field.sizeZ() - 1;
                    Field(i, j, k0 - k) = Field(i, j, k0)
                        +(Field(i, j, k0) - Field(i, j, (k0-1)%Field.sizeZ()))*(-k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }
    // Vectors
    template< class T>
    void BoundaryConditions::SetXVector(Storage3D<T, 0> &Field) const
    {
        if(BC0X == Periodic || BCNX == Periodic)
        {
            for(long int i = -Field.Bcells(); i < 0; i++)
            for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
            for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
            {
                Field(i, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                Field(Field.sizeX() - i - 1, j, k) = Field((-i - 1)%Field.sizeX(), j, k);
            }
            return;
        }

        switch (BC0X)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field((-i - 1)%Field.sizeX(), j, k).Xreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field((-i - 1)%Field.sizeX(), j, k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNX)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(Field.sizeX() - i - 1, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k).Xreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(Field.sizeX() - i - 1, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }

    template< class T>
    void BoundaryConditions::SetXVector(Storage3D<T, 1> &Field) const
    {
        if(BC0X == Periodic || BCNX == Periodic)
        {
            for(long int i = -Field.Bcells(); i < 0; i++)
            for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
            for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
            {
                Field(i, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                Field(Field.sizeX() - i - 1, j, k) = Field((-i - 1)%Field.sizeX(), j, k);
            }
            return;
        }

        switch (BC0X)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                for(long int n = 0; n < Field(0, j, k).size(0); n++)
                {
                    Field(i, j, k)({n}) = Field((-i - 1)%Field.sizeX(), j, k)({n}).Xreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field((-i - 1)%Field.sizeX(), j, k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNX)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                for(long int n = 0; n < Field(Field.sizeX() - 1, j, k).size(0); n++)
                {
                    Field(Field.sizeX() - i - 1, j, k)({n}) = Field((Field.sizeX() + i)%Field.sizeX(), j, k)({n}).Xreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(Field.sizeX() - i - 1, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }
    template< class T>
    void BoundaryConditions::SetXVector(Storage3D<T, 2> &Field) const
    {
        if(BC0X == Periodic || BCNX == Periodic)
        {
            for(long int i = -Field.Bcells(); i < 0; i++)
            for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
            for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
            {
                Field(i, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                Field(Field.sizeX() - i - 1, j, k) = Field((-i - 1)%Field.sizeX(), j, k);
            }
            return;
        }

        switch (BC0X)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                for(long int n1 = 0; n1 < Field(0, j, k).size(0); n1++)
                for(long int n2 = 0; n2 < Field(0, j, k).size(1); n2++)
                {
                    Field(i, j, k)({n1, n2}) = Field((-i - 1)%Field.sizeX(), j, k)({n1, n2}).Xreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field((-i - 1)%Field.sizeX(), j, k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNX)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                for(long int n1 = 0; n1 < Field(Field.sizeX() - 1, j, k).size(0); n1++)
                for(long int n2 = 0; n2 < Field(Field.sizeX() - 1, j, k).size(1); n2++)
                {
                    Field(Field.sizeX() - i - 1, j, k)({n1, n2}) = Field((Field.sizeX() + i)%Field.sizeX(), j, k)({n1, n2}).Xreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(Field.sizeX() - i - 1, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }
    template< class T>
    void BoundaryConditions::SetXVector(Storage3D<T, 3> &Field) const
    {
        if(BC0X == Periodic || BCNX == Periodic)
        {
            for(long int i = -Field.Bcells(); i < 0; i++)
            for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
            for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
            {
                Field(i, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                Field(Field.sizeX() - i - 1, j, k) = Field((-i - 1)%Field.sizeX(), j, k);
            }
            return;
        }

        switch (BC0X)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                for(long int n1 = 0; n1 < Field(0, j, k).size(0); n1++)
                for(long int n2 = 0; n2 < Field(0, j, k).size(1); n2++)
                for(long int n3 = 0; n3 < Field(0, j, k).size(2); n3++)
                {
                    Field(i, j, k)({n1, n2, n3}) = Field((-i - 1)%Field.sizeX(), j, k)({n1, n2, n3}).Xreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field((-i - 1)%Field.sizeX(), j, k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNX)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                for(long int n1 = 0; n1 < Field(Field.sizeX() - 1, j, k).size(0); n1++)
                for(long int n2 = 0; n2 < Field(Field.sizeX() - 1, j, k).size(1); n2++)
                for(long int n3 = 0; n3 < Field(Field.sizeX() - 1, j, k).size(2); n3++)
                {
                    Field(Field.sizeX() - i - 1, j, k)({n1, n2, n3}) = Field((Field.sizeX() + i)%Field.sizeX(), j, k)({n1, n2, n3}).Xreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < 0; i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(Field.sizeX() - i - 1, j, k) = Field((Field.sizeX() + i)%Field.sizeX(), j, k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }
    template< class T>
    void BoundaryConditions::SetYVector(Storage3D<T, 0> &Field) const
    {
        if(BC0Y == Periodic || BCNY == Periodic)
        {
            for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
            for(long int j = -Field.Bcells(); j < 0; j++)
            for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
            {
                Field(i, j, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                Field(i, Field.sizeY() -j - 1, k) = Field(i, (-j - 1)%Field.sizeY(), k);
            }
            return;
        }

        switch (BC0Y)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field(i, (-j - 1)%Field.sizeY(), k).Yreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field(i, (-j - 1)%Field.sizeY(), k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNY)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, Field.sizeY() - j - 1, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k).Yreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, Field.sizeY() - j - 1, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }

    template< class T>
    void BoundaryConditions::SetYVector(Storage3D<T, 1> &Field) const
    {
        if(BC0Y == Periodic || BCNY == Periodic)
        {
            for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
            for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
            {
                Field(i, j, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                Field(i, Field.sizeY() - j - 1, k) = Field(i, (-j - 1)%Field.sizeY(), k);
            }
            return;
        }

        switch (BC0Y)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                for(long int n = 0; n < Field(i, 0, k).size(0); n++)
                {
                    Field(i, j, k)({n}) = Field(i, (-j - 1)%Field.sizeY(), k)({n}).Yreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field(i, (-j - 1)%Field.sizeY(), k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNY)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                for(long int n = 0; n < Field(i, Field.sizeY() - 1, k).size(0); n++)
                {
                    Field(i, Field.sizeY() - j - 1, k)({n}) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k)({n}).Yreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, Field.sizeY() - j - 1, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }
    template< class T>
    void BoundaryConditions::SetYVector(Storage3D<T, 2> &Field) const
    {
        if(BC0Y == Periodic || BCNY == Periodic)
        {
            for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
            for(long int j = -Field.Bcells(); j < 0; j++)
            for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
            {
                Field(i, j, k) = Field(i, Field.sizeY() + j, k);
                Field(i, Field.sizeY() - j - 1, k) = Field(i, (-j - 1)%Field.sizeY(), k);
            }
            return;
        }

        switch (BC0Y)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                for(long int n1 = 0; n1 < Field(i, 0, k).size(0); n1++)
                for(long int n2 = 0; n2 < Field(i, 0, k).size(1); n2++)
                {
                    Field(i, j, k)({n1, n2}) = Field(i, (-j - 1)%Field.sizeY(), k)({n1, n2}).Yreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field(i, (-j - 1)%Field.sizeY(), k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNY)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                for(long int n1 = 0; n1 < Field(i, Field.sizeY() - 1, k).size(0); n1++)
                for(long int n2 = 0; n2 < Field(i, Field.sizeY() - 1, k).size(1); n2++)
                {
                    Field(i, Field.sizeY() - j - 1, k)({n1, n2}) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k)({n1, n2}).Yreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, Field.sizeY() - j - 1, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }
    template< class T>
    void BoundaryConditions::SetYVector(Storage3D<T, 3> &Field) const
    {
        if(BC0Y == Periodic || BCNY == Periodic)
        {
            for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
            for(long int j = -Field.Bcells(); j < 0; j++)
            for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
            {
                Field(i, j, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                Field(i, Field.sizeY() - j - 1, k) = Field(i, (-j - 1)%Field.sizeY(), k);
            }
            return;
        }

        switch (BC0Y)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                for(long int n1 = 0; n1 < Field(i, 0, k).size(0); n1++)
                for(long int n2 = 0; n2 < Field(i, 0, k).size(1); n2++)
                for(long int n3 = 0; n3 < Field(i, 0, k).size(2); n3++)
                {
                    Field(i, j, k)({n1, n2, n3}) = Field(i, (-j - 1)%Field.sizeY(), k)({n1, n2, n3}).Yreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, j, k) = Field(i, (-j - 1)%Field.sizeY(), k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNY)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                for(long int n1 = 0; n1 < Field(i, Field.sizeY() - 1, k).size(0); n1++)
                for(long int n2 = 0; n2 < Field(i, Field.sizeY() - 1, k).size(1); n2++)
                for(long int n3 = 0; n3 < Field(i, Field.sizeY() - 1, k).size(2); n3++)
                {
                    Field(i, Field.sizeY() - j - 1, k)({n1, n2, n3}) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k)({n1, n2, n3}).Yreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < 0; j++)
                for(long int k = -Field.Bcells(); k < Field.sizeZ() + Field.Bcells(); k++)
                {
                    Field(i, Field.sizeY() - j - 1, k) = Field(i, (Field.sizeY() + j)%Field.sizeY(), k);
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }
    template< class T >
    void BoundaryConditions::SetZVector(Storage3D<T, 0> &Field) const
    {
        if(BC0Z == Periodic || BCNZ == Periodic)
        {
            for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
            for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
            for(long int k = -Field.Bcells(); k < 0; k++)
            {
                Field(i, j, k) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (-k - 1)%Field.sizeZ());
            }
            return;
        }

        switch (BC0Z)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, (-k - 1)%Field.sizeZ()).Zreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, (-k - 1)%Field.sizeZ());
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNZ)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ()).Zreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }

    template< class T >
    void BoundaryConditions::SetZVector(Storage3D<T, 1> &Field) const
    {
        if(BC0Z == Periodic || BCNZ == Periodic)
        {
            for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
            for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
            for(long int k = -Field.Bcells(); k < 0; k++)
            {
                Field(i, j, k) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (-k - 1)%Field.sizeZ());
            }
            return;
        }

        switch (BC0Z)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                for(long int n = 0; n < Field(i, j, 0).size(0); n++)
                {
                    Field(i, j, k)({n}) = Field(i, j, (-k - 1)%Field.sizeZ())({n}).Zreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, (-k - 1)%Field.sizeZ());
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNZ)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                for(long int n = 0; n < Field(i, j, Field.sizeZ() - 1).size(0); n++)
                {
                    Field(i, j, Field.sizeZ() - k - 1)({n}) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ())({n}).Zreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }

    template< class T >
    void BoundaryConditions::SetZVector(Storage3D<T, 2> &Field) const
    {
        if(BC0Z == Periodic || BCNZ == Periodic)
        {
            for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
            for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
            for(long int k = -Field.Bcells(); k < 0; k++)
            {
                Field(i, j, k) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (-k - 1)%Field.sizeZ());
            }
            return;
        }

        switch (BC0Z)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                for(long int n1 = 0; n1 < Field(i, j, 0).size(0); n1++)
                for(long int n2 = 0; n2 < Field(i, j, 0).size(1); n2++)
                {
                    Field(i, j, k)({n1, n2}) = Field(i, j, (-k - 1)%Field.sizeZ())({n1, n2}).Zreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, (-k - 1)%Field.sizeZ());
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNZ)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                for(long int n1 = 0; n1 < Field(i, j, Field.sizeZ() - 1).size(0); n1++)
                for(long int n2 = 0; n2 < Field(i, j, Field.sizeZ() - 1).size(1); n2++)
                {
                    Field(i, j, Field.sizeZ() - k - 1)({n1, n2}) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ())({n1, n2}).Zreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }

    template< class T >
    void BoundaryConditions::SetZVector(Storage3D<T, 3> &Field) const
    {
        if(BC0Z == Periodic || BCNZ == Periodic)
        {
            for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
            for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
            for(long int k = -Field.Bcells(); k < 0; k++)
            {
                Field(i, j, k) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                Field(i, j, Field.sizeZ() - k - 1) = Field(i, j, (-k - 1)%Field.sizeZ());
            }
            return;
        }

        switch (BC0Z)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                for(long int n1 = 0; n1 < Field(i, j, 0).size(0); n1++)
                for(long int n2 = 0; n2 < Field(i, j, 0).size(1); n2++)
                for(long int n3 = 0; n3 < Field(i, j, 0).size(2); n3++)
                {
                    Field(i, j, k)({n1, n2, n3}) = Field(i, j, (-k - 1)%Field.sizeZ())({n1, n2, n3}).Zreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, k) = Field(i, j, (-k - 1)%Field.sizeZ());
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }

        switch (BCNZ)
        {
            case NoFlux:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                for(long int n1 = 0; n1 < Field(i, j, Field.sizeZ() - 1).size(0); n1++)
                for(long int n2 = 0; n2 < Field(i, j, Field.sizeZ() - 1).size(1); n2++)
                for(long int n3 = 0; n3 < Field(i, j, Field.sizeZ() - 1).size(2); n3++)
                {
                    Field(i, j, Field.sizeZ() - k - 1)({n1, n2, n3}) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ())({n1, n2, n3}).Zreflected();
                }
                break;
            }
            case Free:
            {
                for(long int i = -Field.Bcells(); i < Field.sizeX() + Field.Bcells(); i++)
                for(long int j = -Field.Bcells(); j < Field.sizeY() + Field.Bcells(); j++)
                for(long int k = -Field.Bcells(); k < 0; k++)
                {
                    Field(i, j, Field.sizeZ() - k -1) = Field(i, j, (Field.sizeZ() + k)%Field.sizeZ());
                }
                break;
            }
            case Fixed:
            default:
            {
                break;
            }
        }
    }
}
#endif