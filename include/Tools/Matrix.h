#ifndef MATRIX_H
#define MATRIX_H

#include "Tools/Includes.h"

namespace opensim
{

template <class T>
class Matrix /// Matrix template class. Can handle any type of values
{
 public:
    Matrix()
    {
        Array = nullptr;
    }
    Matrix(const int nx, const int ny)
    {
    	Allocate(nx, ny);
    }
    T& operator()(const int x, const int y)
    {
#ifdef DEBUG
        if(x >= Size_X or y >= Size_Y)
        {
            std::cout << "Error in Matrix<T>::operator(x,y)\n"
                      << "Access beyond the size of the matrix. Requested position("
                      << x << "," << y << "), max allowed indices "
                      << "(" << Size_X-1 << "," << Size_Y-1
                      << ")\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return Array[Size_Y*x + y];
    }
    void Allocate(const int nx, const int ny)
    {
        Size_X = nx;
        Size_Y = ny;
        Array = new T[nx*ny] ();
    }
    void Reallocate(const int nx, const int ny)
    {
        if (Array != nullptr) delete[] Array;
        Size_X = nx;
        Size_Y = ny;
        Array = new T[nx*ny] ();
    }
    void set(const int x, const int y, const T value)
    {
#ifdef DEBUG
        if(x >= Size_X or y >= Size_Y)
        {
            std::cout << "Error in Matrix<T>::set(x,y)\n"
                      << "Access beyond the size of the matrix. Requested position("
                      << x << "," << y << "), max allowed indices "
                      << "(" << Size_X-1 << "," << Size_Y-1
                      << ")\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        Array[Size_Y*x + y] = value;
    }
    void set_to_value(const T value)
    {
        for (int x = 0; x < Size_X; x++)
        for (int y = 0; y < Size_Y; y++)
        {
            Array[Size_Y*x + y] = value;
        }

    }
    void set_to_zero()
    {
        for (int x = 0; x < Size_X; x++)
        for (int y = 0; y < Size_Y; y++)
        {
            Array[Size_Y*x + y] = T();
        }

    }
    void add(const int x, const int y, const T value)
    {
#ifdef DEBUG
        if(x >= Size_X or y >= Size_Y)
        {
            std::cout << "Error in Matrix<T>::add(x,y)\n"
                      << "Access beyond the size of the matrix. Requested position("
                      << x << "," << y << "), max allowed indices "
                      << "(" << Size_X-1 << "," << Size_Y-1
                      << ")\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        Array[Size_Y*x + y] += value;
    }
    T& get(const int x, const int y) const
    {
#ifdef DEBUG
        if(x >= Size_X or y >= Size_Y)
        {
            std::cout << "Error in Matrix<T>::get(x,y)\n"
                      << "Access beyond the size of the matrix. Requested position("
                      << x << "," << y << "), max allowed indices "
                      << "(" << Size_X-1 << "," << Size_Y-1
                      << ")\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return Array[Size_Y*x + y];
    }
    int size_X() const
    {
        return Size_X;
    }
    int size_Y() const
    {
        return Size_Y;
    }
    T get_min() const
    {
        T min = std::numeric_limits<T>::max();
        for (int x = 0; x < Size_X; x++)
        for (int y = 0; y < Size_Y; y++)
        {
            if (Array[Size_Y*x + y] < min) min = Array[Size_Y*x + y];
        }
        return min;
    }
    T get_max() const
    {
        T max = std::numeric_limits<T>::min();
        for (int x = 0; x < Size_X; x++)
        for (int y = 0; y < Size_Y; y++)
        {
            if (Array[Size_Y*x + y] > max) max = Array[Size_Y*x + y];
        }
        return max;
    }
    bool IsNotAllocated() const
    {
        return (Array == nullptr);
    }
    bool IsAllocated() const
    {
        return !(Array == nullptr);
    }
    std::string print(void) const
    {
        std::stringstream out;
        for(int i = 0; i < Size_X; i++)
        for(int j = 0; j < Size_Y; j++)
        {
            out << "||" << std::setprecision(6) << std::right
                        << std::setw(8) << Array[Size_Y*i+j];
            if (j == Size_Y-1)
                out << "||\n";
            else
                out << " ";
        }
        return out.str();
    };

    Matrix<T>& invert(void)
    {
        if (Size_X != Size_Y)
        {
            std::cout << "Matrix: Can Not Compute Inverse Matrix."
                      << "Matrix:\n" << this->print()
                      << "is not quadratic!!!" << std::endl;
            exit(1);
        }
        std::vector<std::vector<T> > Out;
        std::vector<std::vector<T> > Uni;
        std::vector<T> temp;
        temp.resize(Size_X, 0);
        for(int i = 0; i < Size_X; i++)
        {
            Out.push_back(temp);
            Uni.push_back(temp);
        }

        std::vector<int> indxc(Size_X,0);
        std::vector<int> indxr(Size_X,0);
        std::vector<int> ipiv(Size_X,0);
        int icol = 0;
        int irow = 0;
        T pivinv;
        T dum;
        for(int i = 0; i < Size_X; i++)
        {
            for(int j = 0; j < Size_Y; j++)
            {
                if (i == j)
                {
                    Uni[i][j] = 1.0;
                }
                else
                {
                    Uni[i][j] = 0.0;
                }
                Out[i][j] = Array[Size_Y*i+j];
            }
        }

        for(int i = 0; i < Size_X; i++)
        {
            T big = 0.0;
            for(int j = 0; j < Size_X; j++)
            if(ipiv[j] != 1)
            for(int k = 0; k < Size_Y; k++)
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
                    std::cout << "Matrix: Can Not Compute Inverse Matrix."
                              << "Matrix:\n" << this->print()
                              << "is Singular 1!!!" << std::endl;
                    exit(1);
                }
            };
            ++(ipiv[icol]);
            if(irow != icol)
            {
                for (int l = 0; l < Size_X; l++)
                {
                    T temp = Out[irow][l];
                    Out[irow][l] = Out[icol][l];
                    Out[icol][l] = temp;
                };
                for (int l = 0; l < Size_X; l++)
                {
                    T temp = Uni[irow][l];
                    Uni[irow][l] = Uni[icol][l];
                    Uni[icol][l] = temp;
                };
            };
            indxr[i] = irow;
            indxc[i] = icol;

            if (fabs(Out[icol][icol]) <= DBL_EPSILON)
            {
                std::cout << "Matrix: Can Not Compute Inverse Matrix. Matrix:\n"<< this->print() << "is Singular 2!!!" << std::endl;
                exit(2);
            }
            pivinv = 1.0/Out[icol][icol];
            Out[icol][icol] = 1.0;
            for(int l = 0; l < Size_X; l++) Out[icol][l] *= pivinv;
            for(int l = 0; l < Size_X; l++) Uni[icol][l] *= pivinv;
            for(int ll = 0; ll < Size_X; ll++)
            if(ll != icol)
            {
                dum = Out[ll][icol];
                Out[ll][icol] = 0.0;
                for(int l = 0; l < Size_X; l++) Out[ll][l] -= Out[icol][l]*dum;
                for(int l = 0; l < Size_X; l++) Uni[ll][l] -= Uni[icol][l]*dum;
            }
        }
        for(int l = Size_X-1; l >= 0; l--)
        {
            if(indxr[l] != indxc[l])
            for(int k = 0; k < Size_X; k++)
            {
                T temp = Out[k][indxr[l]];
                Out[k][indxr[l]] = Out[k][indxc[l]];
                Out[k][indxc[l]] = temp;
            };
        }

        for(int i = 0; i < Size_X; i++)
        for(int j = 0; j < Size_Y; j++)
        {
             Array[Size_Y*i+j] = Out[i][j];
        }
        return *this;
    };

    Matrix<T> inverted(void)
    {
        if (Size_X != Size_Y)
        {
            std::cout << "Matrix: Can Not Compute Inverse Matrix."
                      << "Matrix:\n" << this->print()
                      << "is not quadratic!!!" << std::endl;
            exit(1);
        }
        std::vector<std::vector<T> > Out;
        std::vector<std::vector<T> > Uni;
        std::vector<T> temp;
        temp.resize(Size_X, 0);
        for(int i = 0; i < Size_X; i++)
        {
            Out.push_back(temp);
            Uni.push_back(temp);
        }

        std::vector<int> indxc(Size_X,0);
        std::vector<int> indxr(Size_X,0);
        std::vector<int> ipiv(Size_X,0);
        int icol = 0;
        int irow = 0;
        T pivinv;
        T dum;
        for(int i = 0; i < Size_X; i++)
        {
            for(int j = 0; j < Size_Y; j++)
            {
                if (i == j)
                {
                    Uni[i][j] = 1.0;
                }
                else
                {
                    Uni[i][j] = 0.0;
                }
                Out[i][j] = Array[Size_Y*i+j];
            }
        }

        for(int i = 0; i < Size_X; i++)
        {
            T big = 0.0;
            for(int j = 0; j < Size_X; j++)
            if(ipiv[j] != 1)
            for(int k = 0; k < Size_Y; k++)
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
                    std::cout << "Matrix: Can Not Compute Inverse Matrix."
                              << "Matrix:\n" << this->print()
                              << "is Singular 1!!!" << std::endl;
                    exit(1);
                }
            };
            ++(ipiv[icol]);
            if(irow != icol)
            {
                for (int l = 0; l < Size_X; l++)
                {
                    T temp = Out[irow][l];
                    Out[irow][l] = Out[icol][l];
                    Out[icol][l] = temp;
                };
                for (int l = 0; l < Size_X; l++)
                {
                    T temp = Uni[irow][l];
                    Uni[irow][l] = Uni[icol][l];
                    Uni[icol][l] = temp;
                };
            };
            indxr[i] = irow;
            indxc[i] = icol;

            if (fabs(Out[icol][icol]) <= DBL_EPSILON)
            {
                std::cout << "Matrix: Can Not Compute Inverse Matrix. Matrix:\n"<< this->print() << "is Singular 2!!!" << std::endl;
                exit(2);
            }
            pivinv = 1.0/Out[icol][icol];
            Out[icol][icol] = 1.0;
            for(int l = 0; l < Size_X; l++) Out[icol][l] *= pivinv;
            for(int l = 0; l < Size_X; l++) Uni[icol][l] *= pivinv;
            for(int ll = 0; ll < Size_X; ll++)
            if(ll != icol)
            {
                dum = Out[ll][icol];
                Out[ll][icol] = 0.0;
                for(int l = 0; l < Size_X; l++) Out[ll][l] -= Out[icol][l]*dum;
                for(int l = 0; l < Size_X; l++) Uni[ll][l] -= Uni[icol][l]*dum;
            }
        }
        for(int l = Size_X-1; l >= 0; l--)
        {
            if(indxr[l] != indxc[l])
            for(int k = 0; k < Size_X; k++)
            {
                T temp = Out[k][indxr[l]];
                Out[k][indxr[l]] = Out[k][indxc[l]];
                Out[k][indxc[l]] = temp;
            };
        }

        Matrix<T> temp2;
        temp2.Allocate(Size_X, Size_Y);
        for(int i = 0; i < Size_X; i++)
        for(int j = 0; j < Size_Y; j++)
        {
             temp2.set(i,j,Out[i][j]);
        }
        return temp2;
    };

    ~Matrix()
    {
        delete[] Array;
    }
    T* data()
    {
        return Array;
    }
 protected:
 private:
    T*  Array;
    int Size_X;
    int Size_Y;
};
}// namespace opensim
#endif
