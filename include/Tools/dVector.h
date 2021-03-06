#ifndef DVECTOR_H
#define DVECTOR_H

namespace opensim
{

template<int N>
class dVector
{
 public:

    dVector(){};

    dVector(std::initializer_list<double> vecinit)
    {
        for (int i = 0; i < N; i++)
        storage[i] = 0.0;

#ifdef DEBUG
        if (vecinit.size() != N)
        {
            std::cout << "Error in dVector::constructor()\n"
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
        if(i > N-1)
        {
            std::cout << "Error in dVector<int N>::operator[]\n"
                      << "Access beyond storage range. i = "
                      << i << " > N-1"
                      << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i];
    };
    double const& operator[](const int i) const
    {
#ifdef DEBUG
        if(i > N-1)
        {
            std::cout << "Error in dVector<int N>::operator[]\n"
                      << "Access beyond storage range. i = "
                      << i << " > N-1"
                      << "\nTerminating!!!" << std::endl;
            exit(13);
        }
#endif
        return storage[i];
    };
    void set_to_zero(void)
    {
        memset(storage, 0, N*sizeof(double));
    };
    void set_to_value(const double value)
    {
        for(int i = 0; i < N; i++) storage[i] = value;
    };
    dVector<N> operator*(const double m) const
    {
        dVector<N> tmp;
        for (int i = 0; i < N; i++)
        {
            tmp[i] = storage[i]*m;
        }
        return tmp;
    };
    dVector<N>& operator*=(const double m)
    {
        for (int i = 0; i < N; i++)
        {
            storage[i] = storage[i]*m;
        }
        return *this;
    };
    dVector<N> operator+(const double rhs) const
    {
        dVector<N> tmp;
        for (int i = 0; i < N; i++)
        {
            tmp[i] = storage[i] + rhs;
        }
        return tmp;
    };
    dVector<N> operator+(const dVector<N>& rhs) const
    {
        dVector<N> tmp;
        for (int i = 0; i < N; i++)
        {
            tmp[i] = storage[i] + rhs[i];
        }
        return tmp;
    };
    dVector<N>& operator+=(const dVector<N>& rhs)
    {
        for (int i = 0; i < N; i++)
        {
            storage[i] = storage[i] + rhs[i];
        }
        return *this;
    };
    dVector<N> operator-(const double rhs) const
    {
        dVector<N> tmp;
        for (int i = 0; i < N; i++)
        {
            tmp[i] = storage[i] - rhs;
        }
        return tmp;
    };
    dVector<N> operator-(const dVector<N>& rhs) const
    {
        dVector<N> tmp;
        for (int i = 0; i < N; i++)
        {
            tmp[i] = storage[i] - rhs[i];
        }
        return tmp;
    };
    dVector<N>& operator-=(const dVector<N>& rhs)
    {
        for (int i = 0; i < N; i++)
        {
            storage[i] = storage[i] - rhs[i];
        }
        return *this;
    };
    dVector<N>& operator=(const dVector<N>& rhs)
    {
        memmove(data(), rhs.const_data(), N*sizeof(double));
        return *this;
    };
    double sum_of_entries(void) const
    {
        double sum = 0.0;
        for (int i = 0; i < N; i++)
        {
            sum += storage[i];
        }
        return sum;
    };
    double average_entries(void) const
    {
        double sum = 0.0;
        for (int i = 0; i < N; i++)
        {
            sum += storage[i];
        }
        return sum/double(N);
    };
    double* data(void)
    {
        return storage;
    };
    const double* const_data(void) const
    {
        return storage;
    };
    std::string print(void) const
    {
        std::stringstream out;
        out << "||";
        for(int i = 0; i < N-1; i++)
        {
            out << storage[i] << ", ";
        }
        out << storage[N-1] << "||";
        return out.str();
    };

 protected:
 private:
    double storage[N];
};

} // namespace opensim
#endif
