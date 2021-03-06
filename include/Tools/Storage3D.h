#ifndef STORAGE3D_H
#define STORAGE3D_H

#include "Tools/Includes.h"

namespace opensim
{

template <typename A, typename T, typename = void>
struct ClearCaller;

template <typename T>
class has_clear
{
    typedef char one;
    typedef long two;

    template <typename C> static one test( decltype(&C::clear) ) ;
    template <typename C> static two test(...);

public:
    enum { value = sizeof(test<T>(0)) == sizeof(char) };
};

template <typename A, typename T>
class ClearCaller<A, T, typename std::enable_if<has_clear<T>::value>::type>
{
 public:
    static T call(A& self)
    {
        for (long int i = -self.Bcells(); i < self.sizeX() + self.Bcells(); ++i)
        for (long int j = -self.Bcells(); j < self.sizeY() + self.Bcells(); ++j)
        for (long int k = -self.Bcells(); k < self.sizeZ() + self.Bcells(); ++k)
        {
            self(i,j,k).clear();
        }
        return T();
    }
};

template <typename A, typename T>
class ClearCaller<A, T, typename std::enable_if<!has_clear<T>::value>::type> {
public:
    static T call(A& self) {
        std::cerr << "ERROR: Storage3D.h: ClearCaller::Clear() called for non-POD with no clear() method." << std::endl;
        exit(1);
        return T();
    }
};

template <class T, long int Rank>
class Storage3D                                                                 /// 3D storage template class of vector values. Can handle any type of values
{
 public:
    friend class ClearCaller< Storage3D<T, Rank> , T>;

    Storage3D()
    {
        locTensors = nullptr;
        g_level = 1;
    }

    Storage3D<T,Rank> (const long int nx, const long int ny, const long int nz,
            const std::initializer_list<long int> nn, const long int bc, const int grid_level = 1)
    {
        if(grid_level < 1)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " << Rank << ">::Sorage3D()\n"
                      << "Grid level less than 1 is not permitted!" << std::endl;
            std::cerr << "Terminating!!!" << std::endl;
            exit(13);
        }

        std::copy(nn.begin(), nn.end(), Dimensions.begin());

        g_level = grid_level;

        Size_X = nx*g_level;
        Size_Y = ny*g_level;
        Size_Z = nz*g_level;

        b_cells = bc*g_level;

        Size_X_BC = Size_X + 2*b_cells;
        Size_Y_BC = Size_Y + 2*b_cells;
        Size_Z_BC = Size_Z + 2*b_cells;

        locTensors = new Tensor<T, Rank>[Size_X_BC*Size_Y_BC*Size_Z_BC] ();

        for(long int i = 0; i < Size_X_BC*Size_Y_BC*Size_Z_BC; i++)
        {
            locTensors[i].Allocate(Dimensions);
        }
    }

    Storage3D<T,Rank> (const Storage3D<T,Rank>& Field)
    {
        if (this != &Field)
        {
            if (Field.IsAllocated())
            {
                Dimensions = Field.Dimensions;
                g_level = Field.g_level;

                Size_X = Field.Size_X;
                Size_Y = Field.Size_Y;
                Size_Z = Field.Size_Z;
                b_cells = Field.b_cells;
                Size_X_BC = Size_X + 2*b_cells;
                Size_Y_BC = Size_Y + 2*b_cells;
                Size_Z_BC = Size_Z + 2*b_cells;

                locTensors = new Tensor<T, Rank>[Size_X_BC*Size_Y_BC*Size_Z_BC] ();

                for(long int i = 0; i < Size_X_BC*Size_Y_BC*Size_Z_BC; i++)
                {
                    locTensors[i].Allocate(Dimensions);
                }
                const int totsize = Field(0,0,0).size();
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,0,)
                    for (int n = 0; n < totsize; ++n)
                    {
                        locTensors[Index(i,j,k)][n] = Field(i,j,k)[n];
                    }
                OMP_PARALLEL_STORAGE_LOOP_END
            }
            else
            {
                locTensors = nullptr;
                g_level = 1;
            }
        }
    }

    Tensor<T, Rank>& operator()(const long int x, const long int y, const long int z)
    {
#ifdef DEBUG
        if(x > Size_X + b_cells - 1 or y > Size_Y + b_cells - 1 or z > Size_Z + b_cells - 1 or
           x < -b_cells or y < -b_cells or z < -b_cells)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " << Rank << ">::operator()\n"
                 << "Access beyond storage range -> ("
                 << x << "," << y << "," << z << ")" << " is outside storage bounds ["
                 << -b_cells << ", " << -b_cells << ", " << -b_cells << "] and ("
                 << Size_X + b_cells << ", " << Size_Y + b_cells << ", " << Size_Z + b_cells << ")"
                 << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return locTensors[Index(x,y,z)];
    }
    Tensor<T, Rank> const& operator()(const long int x, const long int y,
            const long int z) const
    {
#ifdef DEBUG
        if(x > Size_X + b_cells - 1 or y > Size_Y + b_cells - 1 or z > Size_Z + b_cells - 1 or
           x < -b_cells or y < -b_cells or z < -b_cells)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " << Rank << ">::operator()\n"
                 << "Access beyond storage range -> ("
                 << x << "," << y << "," << z << ")" << " is outside storage bounds ["
                 << -b_cells << ", " << -b_cells << ", " << -b_cells << "] and ("
                 << Size_X + b_cells << ", " << Size_Y + b_cells << ", " << Size_Z + b_cells << ")"
                 << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return locTensors[Index(x,y,z)];
    }

    Tensor<T, Rank> at(const double x, const double y, const double z) const
    {
#ifdef DEBUG
        if(x > Size_X + b_cells - 1 or y > Size_Y + b_cells - 1 or z > Size_Z + b_cells - 1 or
           x < -b_cells or y < -b_cells or z < -b_cells)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " << Rank << ">::operator()\n"
                 << "Access beyond storage range -> ("
                 << x << "," << y << "," << z << ")" << " is outside storage bounds ["
                 << -b_cells << ", " << -b_cells << ", " << -b_cells << "] and ("
                 << Size_X + b_cells << ", " << Size_Y + b_cells << ", " << Size_Z + b_cells << ")"
                 << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif

        long int x0 = floor(x);
        long int y0 = floor(y);
        long int z0 = floor(z);
        double dx = x - x0;
        double dy = y - y0;
        double dz = z - z0;

        Tensor<T, Rank> tempTensor =
              locTensors[Index(x0,y0,z0)]*((1.0 - dx)*(1.0 - dy)*(1.0 - dz)) +
              locTensors[Index(x0+1,y0,z0)]*(dx*(1.0 - dy)*(1.0 - dz)) +
              locTensors[Index(x0,y0+1,z0)]*((1.0 - dx)*dy*(1.0 - dz)) +
              locTensors[Index(x0,y0,z0+1)]*((1.0 - dx)*(1.0 - dy)*dz) +
              locTensors[Index(x0+1,y0+1,z0)]*(dx*dy*(1.0 - dz)) +
              locTensors[Index(x0+1,y0,z0+1)]*(dx*(1.0 - dy)*dz) +
              locTensors[Index(x0,y0+1,z0+1)]*((1.0 - dx)*dy*dz) +
              locTensors[Index(x0+1,y0+1,z0+1)]*(dx*dy*dz);

        return tempTensor;
    }
    Tensor<T, Rank> operator()(const long int x, const long int y, const long int z, const int grid_level)
    {
        switch (grid_level - g_level)
        {
            case 0:
            {
#ifdef DEBUG
                if(x > Size_X + b_cells - 1 or y > Size_Y + b_cells - 1 or z > Size_Z + b_cells - 1 or
                   x < -b_cells or y < -b_cells or z < -b_cells)
                {
                    std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " << Rank << ">::operator()\n"
                         << "Access beyond storage range -> ("
                         << x << "," << y << "," << z << ")" << " is outside of storage bounds ["
                         << -b_cells << ", " << -b_cells << ", " << -b_cells << "] and ("
                         << Size_X + b_cells << ", " << Size_Y + b_cells << ", " << Size_Z + b_cells << ")"
                         << "\nTerminating!!!" << std::endl;
                    exit(13);
                }
#endif
                return locTensors[(Size_Y_BC*(x + b_cells) + y + b_cells)*Size_Z_BC + z + b_cells];
            }
            case -1:
            {
                double i = (x + 0.5)*double(g_level)/double(grid_level) - 0.5;
                double j = (y + 0.5)*double(g_level)/double(grid_level) - 0.5;
                double k = (z + 0.5)*double(g_level)/double(grid_level) - 0.5;
#ifdef DEBUG
                if(i > Size_X + b_cells - 1 or j > Size_Y + b_cells - 1 or k > Size_Z + b_cells - 1 or
                   i < -b_cells or j < -b_cells or k < -b_cells)
                {
                    std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " << Rank << ">::operator()\n"
                         << "Access beyond storage range -> ("
                         << i << "," << j << "," << k << ")" << " is outside storage bounds ["
                         << -b_cells << ", " << -b_cells << ", " << -b_cells << "] and ("
                         << Size_X + b_cells << ", " << Size_Y + b_cells << ", " << Size_Z + b_cells << ")"
                         << "\nTerminating!!!" << std::endl;
                    exit(13);
                }
#endif
                long int x0 = floor(i);
                long int y0 = floor(j);
                long int z0 = floor(k);
                long int x1 = ceil(i);
                long int y1 = ceil(j);
                long int z1 = ceil(k);

                Tensor<T, Rank> tempTensor =
                      (locTensors[Index(x0,y0,z0)] +
                       locTensors[Index(x1,y0,z0)] +
                       locTensors[Index(x0,y1,z0)] +
                       locTensors[Index(x0,y0,z1)] +
                       locTensors[Index(x1,y1,z0)] +
                       locTensors[Index(x1,y0,z1)] +
                       locTensors[Index(x0,y1,z1)] +
                       locTensors[Index(x1,y1,z1)])*0.125;

                return tempTensor;
            }
            case 1:
            {
                double i = (x + 0.5)*double(g_level)/double(grid_level) - 0.5;
                double j = (y + 0.5)*double(g_level)/double(grid_level) - 0.5;
                double k = (z + 0.5)*double(g_level)/double(grid_level) - 0.5;
#ifdef DEBUG
                if(i > Size_X + b_cells - 1 or j > Size_Y + b_cells - 1 or k > Size_Z + b_cells - 1 or
                   i < -b_cells or j < -b_cells or k < -b_cells)
                {
                    std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " << Rank << ">::operator()\n"
                         << "Access beyond storage range -> ("
                         << i << "," << j << "," << k << ")" << " is outside storage bounds ["
                         << -b_cells << ", " << -b_cells << ", " << -b_cells << "] and ("
                         << Size_X + b_cells << ", " << Size_Y + b_cells << ", " << Size_Z + b_cells << ")"
                         << "\nTerminating!!!" << std::endl;
                    exit(13);
                }
#endif
                long int x0 = floor(i);
                long int y0 = floor(j);
                long int z0 = floor(k);
                double dx = i - x0;
                double dy = j - y0;
                double dz = k - z0;

                Tensor<T, Rank> tempTensor =
                      locTensors[Index(x0,y0,z0)]*((1.0 - dx)*(1.0 - dy)*(1.0 - dz)) +
                      locTensors[Index(x0+1,y0,z0)]*(dx*(1.0- dy)*(1.0 - dz)) +
                      locTensors[Index(x0,y0+1,z0)]*((1.0 - dx)*dy*(1.0 - dz)) +
                      locTensors[Index(x0,y0,z0+1)]*((1.0 - dx)*(1.0 - dy)*dz) +
                      locTensors[Index(x0+1,y0+1,z0)]*(dx*dy*(1.0 - dz)) +
                      locTensors[Index(x0+1,y0,z0+1)]*(dx*(1.0 - dy)*dz) +
                      locTensors[Index(x0,y0+1,z0+1)]*((1.0 - dx)*dy*dz) +
                      locTensors[Index(x0+1,y0+1,z0+1)]*(dx*dy*dz);

                return tempTensor;
            }
            default:
            {
                std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " << Rank << ">::operator()\n"
                          << "Wrong grid level requested: <" << grid_level << "> !\n"
                          << "Allowed grid levels are <1> and <2> only!\n"
                          << "Terminating!!!" << std::endl;
                exit(13);
            }
        }
    }

    Tensor<T, Rank> at_grid_level(const double X, const double Y, const double Z, const int grid_level) const
    {
        double x = (X + 0.5)*double(g_level)/double(grid_level) - 0.5;
        double y = (Y + 0.5)*double(g_level)/double(grid_level) - 0.5;
        double z = (Z + 0.5)*double(g_level)/double(grid_level) - 0.5;
#ifdef DEBUG
        if(x > Size_X + b_cells - 1 or y > Size_Y + b_cells - 1 or z > Size_Z + b_cells - 1 or
           x < -b_cells or y < -b_cells or z < -b_cells)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " << Rank << ">::operator()\n"
                 << "Access beyond storage range -> ("
                 << x << "," << y << "," << z << ")" << " is outside storage bounds ["
                 << -b_cells << ", " << -b_cells << ", " << -b_cells << "] and ("
                 << Size_X + b_cells << ", " << Size_Y + b_cells << ", " << Size_Z + b_cells << ")"
                 << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        if(grid_level >= g_level)
        {
            long int x0 = floor(x);
            long int y0 = floor(y);
            long int z0 = floor(z);
            double dx = x - x0;
            double dy = y - y0;
            double dz = z - z0;

            Tensor<T, Rank> tempTensor =
                  locTensors[Index(x0,y0,z0)]*((1.0 - dx)*(1.0 - dy)*(1.0 - dz)) +
                  locTensors[Index(x0+1,y0,z0)]*(dx*(1.0 - dy)*(1.0 - dz)) +
                  locTensors[Index(x0,y0+1,z0)]*((1.0 - dx)*dy*(1.0 - dz)) +
                  locTensors[Index(x0,y0,z0+1)]*((1.0 - dx)*(1.0 - dy)*dz) +
                  locTensors[Index(x0+1,y0+1,z0)]*(dx*dy*(1.0 - dz)) +
                  locTensors[Index(x0+1,y0,z0+1)]*(dx*(1.0 - dy)*dz) +
                  locTensors[Index(x0,y0+1,z0+1)]*((1.0 - dx)*dy*dz) +
                  locTensors[Index(x0+1,y0+1,z0+1)]*(dx*dy*dz);

            return tempTensor;
        }
        else
        {
            long int x0 = floor(x);
            long int y0 = floor(y);
            long int z0 = floor(z);
            long int x1 = ceil(x);
            long int y1 = ceil(y);
            long int z1 = ceil(z);

            Tensor<T, Rank> tempTensor =
                  (locTensors[Index(x0,y0,z0)] +
                   locTensors[Index(x1,y0,z0)] +
                   locTensors[Index(x0,y1,z0)] +
                   locTensors[Index(x0,y0,z1)] +
                   locTensors[Index(x1,y1,z0)] +
                   locTensors[Index(x1,y0,z1)] +
                   locTensors[Index(x0,y1,z1)] +
                   locTensors[Index(x1,y1,z1)])*0.125;

            return tempTensor;
        }
    }

    Tensor<T, Rank>& operator[](const long int idx)
    {
#ifdef DEBUG
         if(idx > Size_X_BC*Size_Y_BC*Size_Z_BC or idx < 0)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " << Rank << ">::operator[]\n"
                      << "Access outside of storage range!\n"
                      << "idx = " << idx << " is outside of allowed range [0, " << Size_X_BC*Size_Y_BC*Size_Z_BC << ")\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return locTensors[idx];
    }
    Tensor<T, Rank>const& operator[](const long int idx) const
    {
#ifdef DEBUG
         if(idx > Size_X_BC*Size_Y_BC*Size_Z_BC or idx < 0)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " << Rank << ">::operator[]\n"
                      << "Access outside of storage range!\n"
                      << "idx = " << idx << " is outside of allowed range [0, " << Size_X_BC*Size_Y_BC*Size_Z_BC << ")\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return locTensors[idx];
    }
    void Allocate(const long int nx, const long int ny, const long int nz,
            const std::initializer_list<long int> nn, const long int bc, const int grid_level = 1)
    {
        if(grid_level < 1)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " <<  Rank << ">::Allocate()\n"
                      << "Grid level less than 1 is not permitted!" << std::endl;
            std::cerr << "Terminating!!!" << std::endl;
            exit(13);
        }

        if(locTensors != nullptr)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " <<  Rank << ">::Allocate()\n"
                      << "Attempt of allocating of a non-empty storage!\n"
                      << "If it is intended, use Reallocate() method instead!" << std::endl;
            std::cerr << "Terminating!!!" << std::endl;
            exit(13);
        }

        std::copy(nn.begin(), nn.end(), Dimensions.begin());

        g_level = grid_level;

        Size_X = nx*g_level;
        Size_Y = ny*g_level;
        Size_Z = nz*g_level;

        b_cells = bc*g_level;

        Size_X_BC = Size_X + 2*b_cells;
        Size_Y_BC = Size_Y + 2*b_cells;
        Size_Z_BC = Size_Z + 2*b_cells;

        locTensors = new Tensor<T, Rank>[Size_X_BC*Size_Y_BC*Size_Z_BC] ();

        for(long int i = 0; i < Size_X_BC*Size_Y_BC*Size_Z_BC; i++)
        {
            locTensors[i].Allocate(Dimensions);
        }
    }
    void Allocate(const Storage3D<T,Rank>& Field)
    {
        if(locTensors != nullptr)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " <<  Rank << ">::Allocate()\n"
                      << "Attempt of allocating of a nonempty storage!\n"
                      << "If it is intended, use Reallocate() method instead!" << std::endl;
            std::cerr << "Terminating!!!" << std::endl;
            exit(13);
        }
        if (this != &Field)
        {
            if (Field.IsAllocated())
            {
                Dimensions = Field.Dimensions;
                g_level = Field.g_level;

                Size_X = Field.Size_X;
                Size_Y = Field.Size_Y;
                Size_Z = Field.Size_Z;
                b_cells = Field.b_cells;
                Size_X_BC = Size_X + 2*b_cells;
                Size_Y_BC = Size_Y + 2*b_cells;
                Size_Z_BC = Size_Z + 2*b_cells;

                locTensors = new Tensor<T, Rank>[Size_X_BC*Size_Y_BC*Size_Z_BC] ();

                for(long int i = 0; i < Size_X_BC*Size_Y_BC*Size_Z_BC; i++)
                {
                    locTensors[i].Allocate(Dimensions);
                }
            }
        }
    }
    void AllocateCopy(const Storage3D<T,Rank>& Field)
    {
        if(locTensors != nullptr)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " <<  Rank << ">::AllocateCopy()\n"
                      << "Attempt of copying to a nonempty storage!\n"
                      << "If it is intended, use assignement operator instead!" << std::endl;
            std::cerr << "Terminating!!!" << std::endl;
            exit(13);
        }

        if (this != &Field)
        {
            if (Field.IsAllocated())
            {
                Dimensions = Field.Dimensions;
                g_level = Field.g_level;

                Size_X = Field.Size_X;
                Size_Y = Field.Size_Y;
                Size_Z = Field.Size_Z;
                b_cells = Field.b_cells;
                Size_X_BC = Size_X + 2*b_cells;
                Size_Y_BC = Size_Y + 2*b_cells;
                Size_Z_BC = Size_Z + 2*b_cells;

                locTensors = new Tensor<T, Rank>[Size_X_BC*Size_Y_BC*Size_Z_BC] ();

                for(long int i = 0; i < Size_X_BC*Size_Y_BC*Size_Z_BC; i++)
                {
                    locTensors[i].Allocate(Dimensions);
                }
                const int totsize = Field(0,0,0).size();
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,0,)
                    for (int n = 0; n < totsize; ++n)
                    {
                        locTensors[Index(i,j,k)][n] = Field(i,j,k)[n];
                    }
                OMP_PARALLEL_STORAGE_LOOP_END
            }
        }
    }
    Storage3D<T,Rank>& operator=(const Storage3D<T,Rank>& Field)
    {
        if (this != &Field)
        {
            if (this->IsAllocated())
            {
                const int totsize = Field(0,0,0).size();
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,0,)
                    for (int n = 0; n < totsize; ++n)
                    {
                        locTensors[Index(i,j,k)][n] = Field(i,j,k)[n];
                    }
                OMP_PARALLEL_STORAGE_LOOP_END
            }
            else
            {
                this->AllocateCopy(Field);
            }
        }
        return *this;
    }
    void Reallocate(const long int nx, const long int ny, const long int nz)
    {

        delete[] locTensors;

        Size_X = nx*g_level;
        Size_Y = ny*g_level;
        Size_Z = nz*g_level;

        Size_X_BC = Size_X + 2*b_cells;
        Size_Y_BC = Size_Y + 2*b_cells;
        Size_Z_BC = Size_Z + 2*b_cells;

        locTensors = new Tensor<T, Rank>[Size_X_BC*Size_Y_BC*Size_Z_BC] ();

        for(long int i = 0; i < Size_X_BC*Size_Y_BC*Size_Z_BC; i++)
        {
            locTensors[i].Allocate(Dimensions);
        }
    }
    bool IsNotAllocated() const
    {
        return (locTensors == nullptr);
    }
    bool IsAllocated() const
    {
        return !(locTensors == nullptr);
    }
    bool IsSize(const long int Nx, const long int Ny, const long int Nz)
    {
        return (Size_X == Nx and Size_Y == Ny and Size_Z == Nz);
    }
    void Remesh(const long int nX, const long int nY, const long int nZ)
    {
        int nx = nX*g_level;
        int ny = nY*g_level;
        int nz = nZ*g_level;

        Tensor<T, Rank>* tempTensors = new Tensor<T, Rank>[(nx + 2*b_cells)*(ny + 2*b_cells)*(nz + 2*b_cells)] ();
        for(long int i = 0; i < (nx + 2*b_cells)*(ny + 2*b_cells)*(nz + 2*b_cells); i++)
        {
            tempTensors[i].Allocate(Dimensions);
        }
        double Xscale = double(Size_X)/double(nx);
        double Yscale = double(Size_Y)/double(ny);
        double Zscale = double(Size_Z)/double(nz);
        #pragma omp parallel for collapse(OMP_COLLAPSE_LOOPS) schedule(dynamic, OMP_DYNAMIC_CHUNKSIZE)
        for(long int x = 0; x < nx; x++)
        for(long int y = 0; y < ny; y++)
        for(long int z = 0; z < nz; z++)
        {
            long int x0 = floor((x - nx*0.5)*Xscale + Size_X * 0.5);
            long int y0 = floor((y - ny*0.5)*Yscale + Size_Y * 0.5);
            long int z0 = floor((z - nz*0.5)*Zscale + Size_Z * 0.5);
            double dx = x*Xscale - x0;
            double dy = y*Yscale - y0;
            double dz = z*Zscale - z0;

            tempTensors[(((ny + 2*b_cells)*(x + b_cells) + y + b_cells)*(nz + 2*b_cells) + z + b_cells)] =
                  locTensors[Index(x0,y0,z0)]*((1.0 - dx)*(1.0 - dy)*(1.0 - dz)) +
                  locTensors[Index(x0+1,y0,z0)]*(dx*(1.0- dy)*(1.0 - dz)) +
                  locTensors[Index(x0,y0+1,z0)]*((1.0 - dx)*dy*(1.0 - dz)) +
                  locTensors[Index(x0,y0,z0+1)]*((1.0 - dx)*(1.0 - dy)*dz) +
                  locTensors[Index(x0+1,y0+1,z0)]*(dx*dy*(1.0 - dz)) +
                  locTensors[Index(x0+1,y0,z0+1)]*(dx*(1.0 - dy)*dz) +
                  locTensors[Index(x0,y0+1,z0+1)]*((1.0 - dx)*dy*dz) +
                  locTensors[Index(x0+1,y0+1,z0+1)]*(dx*dy*dz);
        }
        delete[] locTensors;
        Size_X = nx;
        Size_Y = ny;
        Size_Z = nz;

        Size_X_BC = Size_X + 2*b_cells;
        Size_Y_BC = Size_Y + 2*b_cells;
        Size_Z_BC = Size_Z + 2*b_cells;

        locTensors = tempTensors;
    }
    ~Storage3D()
    {
        delete[] locTensors;
    }

    Tensor<T, Rank>* data(void) const
    {
        return locTensors;
    }

    void Clear(void)
    {
        ClearCaller< Storage3D<T, Rank>, T> K;
        K.call(*this);
    }
    void set_to_value(Tensor<T, Rank>& val)
    {
        if(locTensors != nullptr)
        {
            if(Dimensions == val.Dimensions)
            {
                for(long int idx = 0; idx < Size_X_BC*Size_Y_BC*Size_Z_BC; idx++)
                {
                    locTensors[idx] = val;
                }
            }
            else
            {
                std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " <<  Rank << ">::set_to_value(): Wrong tensor size! \n"
                          << "Terminating!!!" << std::endl;
                exit(1);
            }
        }
        else
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", " <<  Rank << ">::set_to_value() operation on the nonallocated storage! \n"
                      << "Terminating!!!" << std::endl;
            exit(1);
        }
    }
    long int sizeX() const
    {
        return Size_X;
    }
    long int sizeY() const
    {
        return Size_Y;
    }
    long int sizeZ() const
    {
        return Size_Z;
    }
    long int Bcells() const
    {
        return b_cells;
    }

    int GridLevel() const
    {
        return g_level;
    }
    int Dimension() const
    {
        return (Size_X>1) + (Size_Y>1) + (Size_Z>1);
    }
	bool InLimits(const long int x, const long int y, const long int z)
	{
		if (x < -b_cells) return false;
		if (y < -b_cells) return false;
		if (z < -b_cells) return false;
		if (x >= Size_X + b_cells) return false;
		if (y >= Size_Y + b_cells) return false;
		if (z >= Size_Z + b_cells) return false;
		return true;
	}
    void SetNumBcells(const int new_b_cells)
    {
        b_cells = new_b_cells*g_level;
        Reallocate(Size_X, Size_Y, Size_Z);
    }
 protected:
    long int Size_X;
    long int Size_Y;
    long int Size_Z;

    long int Size_X_BC;
    long int Size_Y_BC;
    long int Size_Z_BC;

    long int b_cells;
    int g_level;
    std::array<long int, Rank> Dimensions;
    Tensor<T, Rank>* locTensors;

 private:
    long int Index(const long int x, const long int y, const long int z) const
    {
        return ((Size_Y_BC*(x + b_cells) + y + b_cells)*Size_Z_BC + z + b_cells);
    }
};

template <class T>
class Storage3D<T,0>                                                            /// 3D storage template class specification for Rank == 0. Can handle any type of numerical values
{
 public:
    friend class ClearCaller< Storage3D<T,0>, T>;

    Storage3D<T,0>()
    {
        locData = nullptr;
    }
    Storage3D<T,0> (const Storage3D<T,0>& Field)
    {
        if (this != &Field)
        {
            if (Field.IsAllocated())
            {
                g_level = Field.g_level;

                Size_X = Field.Size_X;
                Size_Y = Field.Size_Y;
                Size_Z = Field.Size_Z;
                b_cells = Field.b_cells;
                Size_X_BC = Size_X + 2*b_cells;
                Size_Y_BC = Size_Y + 2*b_cells;
                Size_Z_BC = Size_Z + 2*b_cells;

                locData = new T[Size_X_BC*Size_Y_BC*Size_Z_BC] ();

                
                OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,Field,b_cells,)
                {
                    locData[Index(i,j,k)] = Field(i,j,k);
                }
                OMP_PARALLEL_STORAGE_LOOP_END
            }
            else
            {
                locData = nullptr;
                g_level = 1;
            }
        }
    }
    Storage3D<T,0> (const long int nx, const long int ny, const long int nz,
            const long int bc, const int grid_level = 1)
    {
        if(grid_level < 1)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::Storage3D(): Grid level less than 1 is not permitted!" << std::endl;
            std::cerr << "Terminating!!!" << std::endl;
            exit(13);
        }

        g_level = grid_level;

        Size_X = nx*g_level;
        Size_Y = ny*g_level;
        Size_Z = nz*g_level;

        b_cells = bc*g_level;

        Size_X_BC = Size_X + 2*b_cells;
        Size_Y_BC = Size_Y + 2*b_cells;
        Size_Z_BC = Size_Z + 2*b_cells;

        locData = new T[Size_X_BC*Size_Y_BC*Size_Z_BC] ();
    }
    T& operator()(const long int x, const long int y, const long int z)
    {
#ifdef DEBUG
        if(x > Size_X + b_cells - 1 or y > Size_Y + b_cells - 1 or z > Size_Z + b_cells - 1 or
           x < -b_cells or y < -b_cells or z < -b_cells)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::operator()\n"
                 << "Access beyond storage range -> ("
                 << x << "," << y << "," << z << ")" << " is outside of storage bounds ["
                 << -b_cells << ", " << -b_cells << ", " << -b_cells << "] and ("
                 << Size_X + b_cells << ", " << Size_Y + b_cells << ", " << Size_Z + b_cells << ")"
                 << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return locData[(Size_Y_BC*(x + b_cells) + y + b_cells)*Size_Z_BC + z + b_cells];
    }
    T const& operator()(const long int x, const long int y, const long int z) const
    {
#ifdef DEBUG
        if(x > Size_X + b_cells - 1 or y > Size_Y + b_cells - 1 or z > Size_Z + b_cells - 1 or
           x < -b_cells or y < -b_cells or z < -b_cells)
        {
            std::cerr << "ERROR: in Storage3D<" << typeid(T).name() << ", 0>::operator()\n"
                 << "Access beyond storage range -> ("
                 << x << "," << y << "," << z << ")" << " is outside of storage bounds ["
                 << -b_cells << ", " << -b_cells << ", " << -b_cells << "] and ("
                 << Size_X + b_cells << ", " << Size_Y + b_cells << ", " << Size_Z + b_cells << ")"
                 << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return locData[(Size_Y_BC*(x + b_cells) + y + b_cells)*Size_Z_BC + z + b_cells];
    }

    T at(const double x, const double y, const double z) const
    {
#ifdef DEBUG
        if(x > Size_X + b_cells - 1 or y > Size_Y + b_cells - 1 or z > Size_Z + b_cells - 1 or
           x < -b_cells or y < -b_cells or z < -b_cells)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::operator()\n"
                 << "Access beyond storage range -> ("
                 << x << "," << y << "," << z << ")" << " is outside of storage bounds ["
                 << -b_cells << ", " << -b_cells << ", " << -b_cells << "] and ("
                 << Size_X + b_cells << ", " << Size_Y + b_cells << ", " << Size_Z + b_cells << ")"
                 << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif

        long int x0 = floor(x);
        long int y0 = floor(y);
        long int z0 = floor(z);
        double dx = x - x0;
        double dy = y - y0;
        double dz = z - z0;

        T tempValue =
              locData[Index(x0,y0,z0)]*((1.0 - dx)*(1.0 - dy)*(1.0 - dz)) +
              locData[Index(x0+1,y0,z0)]*(dx*(1.0- dy)*(1.0 - dz)) +
              locData[Index(x0,y0+1,z0)]*((1.0 - dx)*dy*(1.0 - dz)) +
              locData[Index(x0,y0,z0+1)]*((1.0 - dx)*(1.0 - dy)*dz) +
              locData[Index(x0+1,y0+1,z0)]*(dx*dy*(1.0 - dz)) +
              locData[Index(x0+1,y0,z0+1)]*(dx*(1.0 - dy)*dz) +
              locData[Index(x0,y0+1,z0+1)]*((1.0 - dx)*dy*dz) +
              locData[Index(x0+1,y0+1,z0+1)]*(dx*dy*dz);

        return tempValue;
    }

    T operator()(const long int x, const long int y, const long int z, const int grid_level)
    {
        switch (grid_level - g_level)
        {
            case 0:
            {
#ifdef DEBUG
                if(x > Size_X + b_cells - 1 or y > Size_Y + b_cells - 1 or z > Size_Z + b_cells - 1 or
                   x < -b_cells or y < -b_cells or z < -b_cells)
                {
                    std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::operator()\n"
                         << "Access beyond storage range -> ("
                         << x << "," << y << "," << z << ")" << " is outside of storage bounds ["
                         << -b_cells << ", " << -b_cells << ", " << -b_cells << "] and ("
                         << Size_X + b_cells << ", " << Size_Y + b_cells << ", " << Size_Z + b_cells << ")"
                         << "\nTerminating!!!" << std::endl;
                    exit(13);
                }
#endif
                return locData[(Size_Y_BC*(x + b_cells) + y + b_cells)*Size_Z_BC + z + b_cells];
            }
            case -1:
            {
                double i = (x + 0.5)*double(g_level)/double(grid_level) - 0.5;
                double j = (y + 0.5)*double(g_level)/double(grid_level) - 0.5;
                double k = (z + 0.5)*double(g_level)/double(grid_level) - 0.5;
#ifdef DEBUG
                if(i > Size_X + b_cells - 1 or j > Size_Y + b_cells - 1 or k > Size_Z + b_cells - 1 or
                   i < -b_cells or j < -b_cells or k < -b_cells)
                {
                    std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::operator()\n"
                         << "Access beyond storage range -> ("
                         << i << "," << j << "," << k << ")" << " is outside storage bounds ["
                         << -b_cells << ", " << -b_cells << ", " << -b_cells << "] and ("
                         << Size_X + b_cells << ", " << Size_Y + b_cells << ", " << Size_Z + b_cells << ")"
                         << "\nTerminating!!!" << std::endl;
                    exit(13);
                }
#endif

                long int x0 = floor(i);
                long int y0 = floor(j);
                long int z0 = floor(k);
                long int x1 = ceil(i);
                long int y1 = ceil(j);
                long int z1 = ceil(k);

                T tempValue =
                      (locData[Index(x0,y0,z0)] +
                       locData[Index(x1,y0,z0)] +
                       locData[Index(x0,y1,z0)] +
                       locData[Index(x0,y0,z1)] +
                       locData[Index(x1,y1,z0)] +
                       locData[Index(x1,y0,z1)] +
                       locData[Index(x0,y1,z1)] +
                       locData[Index(x1,y1,z1)])*0.125;

                return tempValue;
            }
            case 1:
            {
                double i = (x + 0.5)*double(g_level)/double(grid_level) - 0.5;
                double j = (y + 0.5)*double(g_level)/double(grid_level) - 0.5;
                double k = (z + 0.5)*double(g_level)/double(grid_level) - 0.5;
#ifdef DEBUG
                if(i > Size_X + b_cells - 1 or j > Size_Y + b_cells - 1 or k > Size_Z + b_cells - 1 or
                   i < -b_cells or j < -b_cells or k < -b_cells)
                {
                    std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::operator()\n"
                         << "Access beyond storage range -> ("
                         << i << "," << j << "," << k << ")" << " is outside storage bounds ["
                         << -b_cells << ", " << -b_cells << ", " << -b_cells << "] and ("
                         << Size_X + b_cells << ", " << Size_Y + b_cells << ", " << Size_Z + b_cells << ")"
                         << "\nTerminating!!!" << std::endl;
                    exit(13);
                }
#endif
                long int x0 = floor(i);
                long int y0 = floor(j);
                long int z0 = floor(k);
                double dx = i - x0;
                double dy = j - y0;
                double dz = k - z0;

                T tempValue =
                      locData[Index(x0,y0,z0)]*((1.0 - dx)*(1.0 - dy)*(1.0 - dz)) +
                      locData[Index(x0+1,y0,z0)]*(dx*(1.0- dy)*(1.0 - dz)) +
                      locData[Index(x0,y0+1,z0)]*((1.0 - dx)*dy*(1.0 - dz)) +
                      locData[Index(x0,y0,z0+1)]*((1.0 - dx)*(1.0 - dy)*dz) +
                      locData[Index(x0+1,y0+1,z0)]*(dx*dy*(1.0 - dz)) +
                      locData[Index(x0+1,y0,z0+1)]*(dx*(1.0 - dy)*dz) +
                      locData[Index(x0,y0+1,z0+1)]*((1.0 - dx)*dy*dz) +
                      locData[Index(x0+1,y0+1,z0+1)]*(dx*dy*dz);

                return tempValue;
            }
            default:
            {
                std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::operator()\n"
                          << "Wrong grid level requested: <" << grid_level << "> !\n"
                          << "Allowed grid levels are <1> and <2> only!\n"
                          << "Terminating!!!" << std::endl;
                exit(13);

            }
        }
        return locData[(Size_Y_BC*(x + b_cells) + y + b_cells)*Size_Z_BC + z + b_cells];
    }
    T at_grid_level(const double X, const double Y, const double Z, const int grid_level) const
    {
        double x = (X + 0.5)*double(g_level)/double(grid_level) - 0.5;
        double y = (Y + 0.5)*double(g_level)/double(grid_level) - 0.5;
        double z = (Z + 0.5)*double(g_level)/double(grid_level) - 0.5;
#ifdef DEBUG
        if(x > Size_X + b_cells - 1 or y > Size_Y + b_cells - 1 or z > Size_Z + b_cells - 1 or
           x < -b_cells or y < -b_cells or z < -b_cells)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::operator()\n"
                 << "Access beyond storage range -> ("
                 << x << "," << y << "," << z << ")" << " is outside of storage bounds ["
                 << -b_cells << ", " << -b_cells << ", " << -b_cells << "] and ("
                 << Size_X + b_cells << ", " << Size_Y + b_cells << ", " << Size_Z + b_cells << ")"
                 << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        if(grid_level >= g_level)
        {
            long int x0 = floor(x);
            long int y0 = floor(y);
            long int z0 = floor(z);
            double dx = x - x0;
            double dy = y - y0;
            double dz = z - z0;

            T tempValue =
                  locData[Index(x0,y0,z0)]*((1.0 - dx)*(1.0 - dy)*(1.0 - dz)) +
                  locData[Index(x0+1,y0,z0)]*(dx*(1.0- dy)*(1.0 - dz)) +
                  locData[Index(x0,y0+1,z0)]*((1.0 - dx)*dy*(1.0 - dz)) +
                  locData[Index(x0,y0,z0+1)]*((1.0 - dx)*(1.0 - dy)*dz) +
                  locData[Index(x0+1,y0+1,z0)]*(dx*dy*(1.0 - dz)) +
                  locData[Index(x0+1,y0,z0+1)]*(dx*(1.0 - dy)*dz) +
                  locData[Index(x0,y0+1,z0+1)]*((1.0 - dx)*dy*dz) +
                  locData[Index(x0+1,y0+1,z0+1)]*(dx*dy*dz);

            return tempValue;
        }
        else
        {
            long int x0 = floor(x);
            long int y0 = floor(y);
            long int z0 = floor(z);
            long int x1 = ceil(x);
            long int y1 = ceil(y);
            long int z1 = ceil(z);

            T tempValue =
                  (locData[Index(x0,y0,z0)] +
                   locData[Index(x1,y0,z0)] +
                   locData[Index(x0,y1,z0)] +
                   locData[Index(x0,y0,z1)] +
                   locData[Index(x1,y1,z0)] +
                   locData[Index(x1,y0,z1)] +
                   locData[Index(x0,y1,z1)] +
                   locData[Index(x1,y1,z1)])*0.125;

            return tempValue;
        }
    }
    T& operator[](long int idx)
    {
#ifdef DEBUG
        if(idx > Size_X_BC*Size_Y_BC*Size_Z_BC - 1 or idx < 0)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::operator[]\n"
                      << "Access outside of storage range!\n"
                      << "idx = " << idx << " is outside of allowed range [0, " << Size_X_BC*Size_Y_BC*Size_Z_BC << ")\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return locData[idx];
    }
    T const& operator[](const long int idx) const
    {
#ifdef DEBUG
        if(idx > Size_X_BC*Size_Y_BC*Size_Z_BC - 1 or idx < 0)
        {
            std::cerr << "ERROR: Access outside storage range in Storage3D<" << typeid(T).name() << ", 0>::operator[]\n"
                      << "idx = " << idx << " is outside of allowed range [0, " << Size_X_BC*Size_Y_BC*Size_Z_BC << ")\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return locData[idx];
    }
    void Allocate(const long int nx, const long int ny, const long int nz,
            const long int bc, const int grid_level = 1)
    {
        if(grid_level < 1)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::Allocate()\n"
                      << "Grid level less than 1 is not permitted!\n"
                      << "Terminating!!!" << std::endl;
            exit(13);
        }
        if(locData != nullptr)
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::Allocate()\n"
                      << "Attempt of reallocation of a non-empty storage!\n"
                      << "If it is intended, use Reallocate() method instead!\n"
                      << "Terminating!!!" << std::endl;
            exit(13);
        }

        g_level = grid_level;

        Size_X = nx*g_level;
        Size_Y = ny*g_level;
        Size_Z = nz*g_level;

        b_cells = bc*g_level;

        Size_X_BC = Size_X + 2*b_cells;
        Size_Y_BC = Size_Y + 2*b_cells;
        Size_Z_BC = Size_Z + 2*b_cells;

        locData = new T[Size_X_BC*Size_Y_BC*Size_Z_BC] ();
    }
    void Reallocate(const long int nx, const long int ny, const long int nz)
    {
        delete[] locData;

        Size_X = nx*g_level;
        Size_Y = ny*g_level;
        Size_Z = nz*g_level;

        Size_X_BC = Size_X + 2*b_cells;
        Size_Y_BC = Size_Y + 2*b_cells;
        Size_Z_BC = Size_Z + 2*b_cells;
        locData = new T[Size_X_BC*Size_Y_BC*Size_Z_BC] ();
    }
    bool IsNotAllocated() const
    {
        return (locData == nullptr);
    }
    bool IsAllocated() const
    {
        return !(locData == nullptr);
    }
    bool IsSize(const long int Nx, const long int Ny, const long int Nz)
    {
        return (Size_X == Nx and Size_Y == Ny and Size_Z == Nz);
    }
    int Dimension() const
    {
        return (Size_X>1) + (Size_Y>1) + (Size_Z>1);
    }
    void Remesh(const long int nX, const long int nY, const long int nZ)
    {
        int nx = nX*g_level;
        int ny = nY*g_level;
        int nz = nZ*g_level;

        T* tempArray = new T[(nx + 2*b_cells)*(ny + 2*b_cells)*(nz + 2*b_cells)] ();

        double Xscale = double(Size_X)/double(nx);
        double Yscale = double(Size_Y)/double(ny);
        double Zscale = double(Size_Z)/double(nz);

        #pragma omp parallel for collapse(OMP_COLLAPSE_LOOPS) schedule(dynamic, OMP_DYNAMIC_CHUNKSIZE)
        for(long int x = 0; x < nx; x++)
        for(long int y = 0; y < ny; y++)
        for(long int z = 0; z < nz; z++)
        {
            long int x0 = floor((x - nx*0.5)*Xscale + Size_X * 0.5);
            long int y0 = floor((y - ny*0.5)*Yscale + Size_Y * 0.5);
            long int z0 = floor((z - nz*0.5)*Zscale + Size_Z * 0.5);
            double dx = x*Xscale - x0;
            double dy = y*Yscale - y0;
            double dz = z*Zscale - z0;

            tempArray[(((ny + 2*b_cells)*(x + b_cells) + y + b_cells)*(nz + 2*b_cells) + z + b_cells)] =
                locData[Index(x0,y0,z0)]*((1.0 - dx)*(1.0 - dy)*(1.0 - dz)) +
                locData[Index(x0+1,y0,z0)]*(dx*(1.0- dy)*(1.0 - dz)) +
                locData[Index(x0,y0+1,z0)]*((1.0 - dx)*dy*(1.0 - dz)) +
                locData[Index(x0,y0,z0+1)]*((1.0 - dx)*(1.0 - dy)*dz) +
                locData[Index(x0+1,y0+1,z0)]*(dx*dy*(1.0 - dz)) +
                locData[Index(x0+1,y0,z0+1)]*(dx*(1.0 - dy)*dz) +
                locData[Index(x0,y0+1,z0+1)]*((1.0 - dx)*dy*dz) +
                locData[Index(x0+1,y0+1,z0+1)]*(dx*dy*dz);
        }
        delete[] locData;
        Size_X = nx;
        Size_Y = ny;
        Size_Z = nz;

        Size_X_BC = Size_X + 2*b_cells;
        Size_Y_BC = Size_Y + 2*b_cells;
        Size_Z_BC = Size_Z + 2*b_cells;

        locData = tempArray;
    }

    bool rotate(const long int newdimx, const long int newdimy, const long int newdimz)
    {
        if ((newdimx == newdimy) or (newdimx == newdimz) or (newdimy == newdimz))
        {
            std::cerr << "flip_storage: Double rotation axis.";
            return false;
        }
        if (newdimx < 0 or newdimx > 2)
        {
            std::cerr << "flip_storage: Wrong parameter newdimx.";
            return false;
        }
        if (newdimy < 0 or newdimy > 2)
        {
            std::cerr << "flip_storage: Wrong parameter newdimy.";
            return false;
        }
        if (newdimz < 0 or newdimz > 2)
        {
            std::cerr << "flip_storage: Wrong parameter newdimz.";
            return false;
        }

        long int newdim[3] = {newdimx, newdimy, newdimz};
        long int oldSize[3] = {Size_X, Size_Y, Size_Z};

        // newdimx == 0; newdimy == 1; newdimz == 2
        if (newdimx == 0 and newdimy == 1 and newdimz == 2) return false;

        T* tempArray = new T[(Size_X + 2*b_cells)*(Size_Y + 2*b_cells)*(Size_Z + 2*b_cells)] ();

        // dimx == 0; dimy == 2; dimz == 1 (rotates around postive x)
        if (newdimx == 0 and newdimy == 2 and newdimz == 1)
        {
            for (int i = 0; i < Size_X; i++)
            for (int j = 0; j < Size_Y; j++)
            for (int k = 0; k < Size_Z; k++)
            {
                tempArray[(Size_Z_BC*(i + b_cells) + Size_Z - 1 - k + b_cells)*Size_Y_BC + j + b_cells]
                          = locData[Index(i,j,k)];
            }
        }
        // dimx == 2; dimy == 1; dimz == 0 (rotates around postive y)
        if (newdimx == 2 and newdimy == 1 and newdimz == 0)
        {
            for (int i = 0; i < Size_X; i++)
            for (int j = 0; j < Size_Y; j++)
            for (int k = 0; k < Size_Z; k++)
            {
                tempArray[(Size_Y_BC*(k + b_cells) + j + b_cells)*Size_X_BC + (Size_X - i - 1) + b_cells]
                          = locData[Index(i,j,k)];
            }
        }
        // dimx == 1; dimy == 0; dimz == 2 (rotates around postive z)
        if (newdimx == 1 and newdimy == 0 and newdimz == 2)
        {
            for (int i = 0; i < Size_X; i++)
            for (int j = 0; j < Size_Y; j++)
            for (int k = 0; k < Size_Z; k++)
            {
                tempArray[(Size_X_BC*(j + b_cells) + i + b_cells)*Size_Z_BC + k + b_cells]
                          = locData[Index(i,j,k)];
            }
        }

        delete[] locData;
        Size_X = oldSize[newdim[0]];
        Size_Y = oldSize[newdim[1]];
        Size_Z = oldSize[newdim[2]];

        Size_X_BC = Size_X + 2*b_cells;
        Size_Y_BC = Size_Y + 2*b_cells;
        Size_Z_BC = Size_Z + 2*b_cells;
        locData = tempArray;

        return true;
    };

    ~Storage3D<T,0>()
    {
        delete[] locData;
    }
    Storage3D<T,0>& operator=(const Storage3D<T,0>& locStorage3D)
    {
        if(locData != nullptr and locStorage3D.IsAllocated())
        {        
            if(std::is_arithmetic<T>::value)
            {
                if (locStorage3D.Size_X_BC == Size_X_BC and
                    locStorage3D.Size_Y_BC == Size_Y_BC and
                    locStorage3D.Size_Z_BC == Size_Z_BC)
                {
                    memcpy(locData, locStorage3D.data(), sizeof(T)*Size_X_BC*Size_Y_BC*Size_Z_BC);
                }
                else
                {
                    std::cerr << "ERROR: Wrong storage size in Storage3D<" << typeid(T).name() << ", 0>::operator=()!\n"
                              << "Terminating!!!" << std::endl;
                    exit(1);
                }
            }
            else
            {
                if (locStorage3D.Size_X_BC == Size_X_BC and
                    locStorage3D.Size_Y_BC == Size_Y_BC and
                    locStorage3D.Size_Z_BC == Size_Z_BC)
                {
                    OMP_PARALLEL_STORAGE_LOOP_BEGIN(i,j,k,locStorage3D,locStorage3D.Bcells(),)
                    {
                        locData[Index(i,j,k)] = locStorage3D(i,j,k);
                    }
                    OMP_PARALLEL_STORAGE_LOOP_END
                }
                else
                {
                    std::cerr << "ERROR: Wrong storage size in Storage3D<" << typeid(T).name() << ", 0>::operator=()!\n"
                              << "Terminating!!!" << std::endl;
                    exit(1);
                }
            }
        }        
        else if (locData == nullptr and locStorage3D.IsAllocated())
        {
            Storage3D<T,0>(locStorage3D);
        }
        return *this;
    }
    void Clear(void)
    {
        if(locData != nullptr and std::is_arithmetic<T>::value)
        {
            memset(locData, 0, sizeof(T)*Size_X_BC*Size_Y_BC*Size_Z_BC);
        }
        else
        {
            ClearCaller< Storage3D<T,0> ,T> K;
            K.call(*this);
        }
    }
    void set_to_value(T val)
    {
        if(locData != nullptr)
        {
            for(long int idx = 0; idx < Size_X_BC*Size_Y_BC*Size_Z_BC; idx++)
            {
                locData[idx] = val;
            }
        }
        else
        {
            std::cerr << "ERROR: Storage3D<" << typeid(T).name() << ", 0>::set_to_value() operation on the nonallocated storage! \n"
                      << "Terminating!!!" << std::endl;
            exit(1);
        }
    }
    T* data(void)
    {
        return locData;
    }
    T* data(void) const
    {
        return locData;
    }
    long int sizeX() const
    {
        return Size_X;
    }
    long int sizeY() const
    {
        return Size_Y;
    }
    long int sizeZ() const
    {
        return Size_Z;
    }
    long int Bcells() const
    {
        return b_cells;
    }
    long int tot_size() const
    {
        return Size_X_BC*Size_Y_BC*Size_Z_BC;
    }
    int GridLevel() const
    {
        return g_level;
    }
    void SetNumBcells(const int new_b_cells)
    {
        b_cells = new_b_cells*g_level;
        Reallocate(Size_X, Size_Y, Size_Z);
    }
    std::string print(void) const
    {
        std::stringstream out;
        for (long int k = 0; k < Size_Z; ++k)
        for (long int j = 0; j < Size_Y; ++j)
        for (long int i = 0; i < Size_X; ++i)
        {
            out << " " << locData[Index(i,j,k)] << std::endl;
        }
        return out.str();
    };
    std::string print_slice(const std::string dim, const long int scoord) const
    {
        std::stringstream out;

        if (dim == "z" or dim == "Z")
        if (scoord <= Size_Z)
        {
            for (long int j = 0; j < Size_Y; ++j)
            {
                for (long int i = 0; i < Size_X; ++i)
                {
                    out << " " << locData[Index(i,j,scoord)];
                }
                out << std::endl;
            }
        }
        if (dim == "y" or dim == "Y")
        if (scoord <= Size_Y)
        {
            for (long int k = 0; k < Size_Z; ++k)
            {
                for (long int i = 0; i < Size_X; ++i)
                {
                    out << " " << locData[Index(i,scoord,k)];
                }
                out << std::endl;
            }
        }
        if (dim == "x" or dim == "X")
        if (scoord <= Size_X)
        {
            for (long int k = 0; k < Size_Z; ++k)
            {
                for (long int j = 0; j < Size_Y; ++j)
                {
                    out << " " << locData[Index(scoord,j,k)];
                }
                out << std::endl;
            }
        }
        return out.str();
    };
    std::string print_subset(const long int xmin, const long int xmax,
                              const long int ymin, const long int ymax,
                              const long int zmin, const long int zmax) const
    {
        std::stringstream out;
        if(xmin < -b_cells or xmin > xmax or xmin > sizeX() + b_cells) {out << "Storage3D.print() - bad index"; return out.str();}
        if(ymin < -b_cells or xmin > xmax or xmin > sizeY() + b_cells) {out << "Storage3D.print() - bad index"; return out.str();}
        if(zmin < -b_cells or zmin > xmax or zmin > sizeZ() + b_cells) {out << "Storage3D.print() - bad index"; return out.str();}
        if(xmax < -b_cells or xmax < xmin or xmax > sizeX() + b_cells) {out << "Storage3D.print() - bad index"; return out.str();}
        if(ymax < -b_cells or ymax < ymin or ymax > sizeY() + b_cells) {out << "Storage3D.print() - bad index"; return out.str();}
        if(zmax < -b_cells or zmax < zmin or zmax > sizeZ() + b_cells) {out << "Storage3D.print() - bad index"; return out.str();}

        for (long int k = zmin; k < zmax; ++k)
        {
            for (long int j = ymin; j < ymax; ++j)
            {
                for (long int i = xmin; i < xmax; ++i)
                {
                    out << " " << locData[Index(i,j,k)];
                }
                out << std::endl;
            }
            out << std::endl;
        }
        return out.str();
    };
    void WriteVTK(const std::string filename, int nbcells = 0) const
    {
        if (nbcells > b_cells) nbcells = b_cells;

        if(std::is_arithmetic<T>::value)
        {
            std::stringstream outbufer;

            outbufer << "# vtk DataFile Version 3.0\n";
            outbufer << "Storage3D\n";
            outbufer << "ASCII\n";
            outbufer << "DATASET STRUCTURED_GRID\n";
            outbufer << "DIMENSIONS " << Size_X + 2*nbcells << " "
                                      << Size_Y + 2*nbcells << " "
                                      << Size_Z + 2*nbcells << "\n";
            outbufer << "POINTS " <<  (Size_Z + 2*nbcells)*(Size_Y + 2*nbcells)
                            *(Size_X+2*nbcells) << " int\n";

            for (int k = -nbcells; k < Size_Z + nbcells; ++k)
            for (int j = -nbcells; j < Size_Y + nbcells; ++j)
            for (int i = -nbcells; i < Size_X + nbcells; ++i)
            {
                outbufer << i << " " << j << " " << k << "\n";
            }
            outbufer << " \n";
            outbufer << "POINT_DATA " << (Size_Z + 2*nbcells)*(Size_Y + 2*nbcells)
                    *(Size_X+2*nbcells) << " \n";
            outbufer << "SCALARS " << filename << " double\n";
            outbufer << "LOOKUP_TABLE default\n";

            for (int k = 0; k < Size_Z_BC; ++k)
            for (int j = 0; j < Size_Y_BC; ++j)
            for (int i = 0; i < Size_X_BC; ++i)
            if (k > nbcells or k < Size_Z_BC - nbcells)
            if (j > nbcells or j < Size_Y_BC - nbcells)
            if (i > nbcells or i < Size_X_BC - nbcells)
            {
                outbufer << locData[Index(i-nbcells,j-nbcells,k-nbcells)] << " ";
            }
            outbufer << std::endl;

            std::string FileName = "VTK/" + filename + ".vtk";

            std::ofstream vtk_file(FileName.c_str());
            vtk_file << outbufer.rdbuf();
            vtk_file.close();
        }
    }
	bool InLimits(const long int x, const long int y, const long int z)
	{
		if (x < -b_cells) return false;
		if (y < -b_cells) return false;
		if (z < -b_cells) return false;
		if (x >= Size_X + b_cells) return false;
		if (y >= Size_Y + b_cells) return false;
		if (z >= Size_Z + b_cells) return false;
		return true;
	}

 protected:
 private:
    long int Size_X;
    long int Size_Y;
    long int Size_Z;

    long int Size_X_BC;
    long int Size_Y_BC;
    long int Size_Z_BC;

    long int b_cells;
    int g_level;
    T* locData;

    long int Index(const long int x, const long int y, const long int z) const
    {
        return ((Size_Y_BC*(x + b_cells) + y + b_cells)*Size_Z_BC + z + b_cells);
    }
};

}// namespace opensim
#endif