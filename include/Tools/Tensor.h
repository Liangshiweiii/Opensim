#ifndef TENSOR_H
#define TENSOR_H

namespace opensim{
    template <class T, long int Rank>
    class Tensor                                                                    /// Tensor template class. Can handle any type of values but POD is preferred
    {
        public:
        Tensor()
        {
            locData = nullptr;
            Dimensions.fill(0);
            allocated = false;
            totSize = 0;
        }
        Tensor(std::array<long int, Rank> locDimensions)
        {
            Dimensions = locDimensions;
            totSize = 1;
            for(long int n = 0; n < Dimensions.size(); n++)
            {
                totSize *= Dimensions[n];
            }
            if(totSize)
            {
                locData = new T[totSize] ();
                allocated = true;
            }
            else
            {
                locData = nullptr;
                allocated = false;
            }
        }
        Tensor(std::initializer_list<long int> locDimensions)
        {
            std::copy(locDimensions.begin(), locDimensions.end(), Dimensions.begin());
            totSize = 1;
            for(unsigned int n = 0; n < Dimensions.size(); n++)
            {
                totSize *= Dimensions[n];
            }

            if(totSize)
            {
                locData = new T[totSize] ();
                allocated = true;
            }
            else
            {
                locData = nullptr;
                allocated = false;
            }
        }
        Tensor(T* data_ptr, std::array<long int, Rank> locDimensions)
        {
            locData = data_ptr;
            allocated = false;
            Dimensions = locDimensions;
            totSize = 1;
            for(unsigned int n = 0; n < Dimensions.size(); n++)
            {
                totSize *= Dimensions[n];
            }
        }
        Tensor(const Tensor<T, Rank>& locTensor)
        {
            Dimensions = locTensor.Dimensions;
            totSize = locTensor.totSize;

            if(totSize)
            {
                locData = new T[totSize] ();
                memmove(locData, locTensor.data(), sizeof(T)*totSize);
                allocated = true;
            }
            else
            {
                locData = nullptr;
                allocated = false;
            }
        }
        T& operator()(std::initializer_list<long int> locPosition)
        {
    #ifdef DEBUG

            if(locPosition.size() != rank())
            {
                std::cout << "Error in Tensor<T," << Rank << ">::operator()\n"
                    << "Access beyond the rank of the tensor. Requested rank = "
                    << locPosition.size() << " > tensor rank = " << Rank
                    << "\nTerminating!!!" << std::endl;
                exit(13);
            }
            else
            {
                std::array<long int, Rank> position;
                std::copy(locPosition.begin(), locPosition.end(), position.begin());
                for(long int n = 0; n < rank(); n++)
                if(position[n] >= Dimensions[n])
                {
                    std::cout << "Error in Tensor<T," << Rank << ">::operator()\n"
                            << "Access beyond the size of the tensor. Requested position[" << n << "] = "
                            << position[n] << " > allowed size of "
                            << Dimensions[n]
                            << "\nTerminating!!!" << std::endl;
                    exit(13);
                }
            }
    #endif
            return locData[Index(locPosition)];
        }
        T const& operator()(const std::initializer_list<long int> locPosition) const
        {
    #ifdef DEBUG

            if(locPosition.size() != rank())
            {
                std::cout << "Error in Tensor<T," << Rank << ">::operator()\n"
                    << "Access beyond the rank of the tensor. Requested rank = "
                    << locPosition.size() << " > tensor rank = " << Rank
                    << "\nTerminating!!!" << std::endl;
                exit(13);
            }
            else
            {
                std::array<long int, Rank> position;
                std::copy(locPosition.begin(), locPosition.end(), position.begin());
                for(long int n = 0; n < rank(); n++)
                if(position[n] >= Dimensions[n])
                {
                    std::cout << "Error in Tensor<T," << Rank << ">::operator()\n"
                            << "Access beyond the size of the tensor. Requested position[" << n << "] = "
                            << position[n] << " >= allowed size of "
                            << Dimensions[n]
                            << "\nTerminating!!!" << std::endl;
                    exit(13);
                }
            }
    #endif
            return locData[Index(locPosition)];
        }
        T& operator[](const unsigned long int position)
        {
    #ifdef DEBUG
            if(position >= totSize)
            {
                std::cout << "Error in Tensor<T," << Rank << ">::operator[]\n"
                        << "Access beyond the total size of the tensor. Requested position = "
                        << position << " >= tensor total size = " << totSize
                        << "\nTerminating!!!" << std::endl;
                exit(13);
            }
    #endif
            return locData[position];
        }
        T const& operator[](const unsigned long int position) const
        {
    #ifdef DEBUG
            if(position >= totSize)
            {
                std::cout << "Error in Tensor<T," << Rank << ">::operator[]\n"
                        << "Access beyond the total size of the tensor. Requested position = "
                        << position << " > tensor total size = " << totSize
                        << "\nTerminating!!!" << std::endl;
                exit(13);
            }
    #endif
            return locData[position];
        }

        void Assign(T* data_ptr, std::array<long int, Rank> locDimensions)
        {
            locData = data_ptr;
            allocated = false;
            Dimensions = locDimensions;
            totSize = 1;
            for(unsigned int n = 0; n < Dimensions.size(); n++)
            {
                totSize *= Dimensions[n];
            }
        }

        void Allocate(std::array<long int, Rank> locDimensions)
        {
            if(!allocated)
            {
                Dimensions = locDimensions;
                totSize = 1;
                for(unsigned int n = 0; n < Dimensions.size(); n++)
                {
                    totSize *= Dimensions[n];
                }
                if(totSize)
                {
                    locData = new T[totSize] ();
                    allocated = true;
                }
                else
                {
                    locData = nullptr;
                    allocated = false;
                }
            }
            else
            {
                std::cout << "ERROR: Tensor<T," << Rank << ">: Allocation attempt of already allocated Tensor! Reallocate() should be used instead!"
                        << std::endl;
                exit(133);
            }
        }

        void Allocate(std::initializer_list<long int> locDimensions)
        {
            if(!allocated)
            {
                std::copy(locDimensions.begin(), locDimensions.end(), Dimensions.begin());
                totSize = 1;
                for(unsigned int n = 0; n < Dimensions.size(); n++)
                {
                    totSize *= Dimensions[n];
                }
                if(totSize)
                {
                    locData = new T[totSize] ();
                    allocated = true;
                }
                else
                {
                    locData = nullptr;
                    allocated = false;
                }
            }
            else
            {
                std::cout << "ERROR: Tensor<T," << Rank << ">: Allocation attempt of already allocated Tensor! Reallocate() should be used instead!"
                        << std::endl;
                exit(133);
            }
        }
        void Reallocate(std::initializer_list<long int> locDimensions)
        {
            if(allocated)
            {
                delete[] locData;
            }

            std::copy(locDimensions.begin(), locDimensions.end(), Dimensions.begin());
            totSize = 1;
            for(unsigned int n = 0; n < Dimensions.size(); n++)
            {
                totSize *= Dimensions[n];
            }

            if(totSize)
            {
                locData = new T[totSize] ();
                allocated = true;
            }
            else
            {
                locData = nullptr;
                allocated = false;
            }
        }
        void set_to_value(T value)
        {
            if(allocated)
            for (unsigned int i = 0; i < totSize; i++)
            {
                locData[i] = value;
            }
        }
        void set_to_zero(void)
        {
            if(allocated)
            for (unsigned int i = 0; i < totSize; i++)
            {
                locData[i] = T();
            }
        }
        unsigned int size(unsigned int n) const
        {
            if(n < Dimensions.size())
            {
                return Dimensions[n];
            }
            else
            {
                return 0;
            }
        }
        unsigned int size(void) const
        {
            return totSize;
        }
        unsigned int rank(void) const
        {
            return Dimensions.size();
        }
        ~Tensor()
        {
            if(allocated)
            {
                delete[] locData;
            }
        }
        const T* data(void) const
        {
            return locData;
        }
        T* data(void)
        {
            return locData;
        }
        Tensor<T, Rank>& operator=(const Tensor<T, Rank>& locTensor)
        {
            if (this == &locTensor)
            {
                return *this;
            }
            if (!locTensor.allocated)
            {
                if(allocated)
                {
                    delete[] locData;
                }
                allocated = false;
                return *this;
            }
            if (!allocated)
            {
                Dimensions = locTensor.Dimensions;
                totSize = 1;
                for(unsigned long int n = 0; n < Dimensions.size(); n++)
                {
                    totSize *= Dimensions[n];
                }
                if(totSize)
                {
                    locData = new T[totSize] ();
                    allocated = true;
                }
            }
            if (locTensor.Dimensions == Dimensions)
            {
                memmove(locData, locTensor.data(), sizeof(T)*totSize);
            }
            else
            {
                std::cout << "ERROR: Tensor<T," << Rank << ">: Different tensor size in assignment operator!"
                        << std::endl;
                exit(133);
            }
                return *this;
        }

        template<typename T2>
        Tensor<T, Rank>& operator+=(const Tensor<T2, Rank>& locTensor)
        {
            if ((not allocated) and (locTensor.allocated))
            {
                Dimensions = locTensor.Dimensions;
                totSize = 1;
                for(unsigned long int n = 0; n < Dimensions.size(); n++)
                {
                    totSize *= Dimensions[n];
                }
                if(totSize)
                {
                    locData = new T[totSize] ();
                    allocated = true;
                }
            }
            if (locTensor.Dimensions == Dimensions)
            {
                for(unsigned long int n = 0; n < totSize; n++)
                {
                    locData[n] += (T)locTensor[n];
                }
            }
            else
            {
                std::cout << "ERROR: Tensor<T," << Rank
                    << ">: Different tensor size in += operator!" << std::endl;

                std::cout << "Dimensions lside: ";
                    for (int i = 0; i < Rank; i++)
                        std::cout << Dimensions[i] << " ";
                std::cout << std::endl;
                std::cout << "Dimensions rside: ";
                    for (int i = 0; i < Rank; i++)
                        std::cout << locTensor.Dimensions[i] << " ";
                std::cout << std::endl;

                exit(133);
            }
            return *this;
        }

        template<typename T2>
        Tensor<T, Rank>& operator-=(const Tensor<T2, Rank>& locTensor)
        {
            if ((not allocated) and (locTensor.allocated))
            {
                Dimensions = locTensor.Dimensions;
                totSize = 1;
                for(unsigned long int n = 0; n < Dimensions.size(); n++)
                {
                    totSize *= Dimensions[n];
                }
                if(totSize)
                {
                    locData = new T[totSize] ();
                    allocated = true;
                }
            }
            if (locTensor.Dimensions == Dimensions)
            {
                for(unsigned long int n = 0; n < totSize; n++)
                {
                    locData[n] += (T)locTensor[n];
                }
            }
            else
            {
                std::cout << "ERROR: Tensor<T," << Rank
                    << ">: Different tensor size in -= operator!"
                    << std::endl;

                std::cout << "Dimensions lside: ";
                    for (int i = 0; i < Rank; i++)
                        std::cout << Dimensions[i] << " ";
                std::cout << std::endl;
                std::cout << "Dimensions rside: ";
                    for (int i = 0; i < Rank; i++)
                        std::cout << locTensor.Dimensions[i] << " ";
                std::cout << std::endl;

                exit(133);
            }
            return *this;
        }

        template<typename T2>
        Tensor<T, Rank>& operator/=(const T2 number)
        {
            if (allocated)
            for(unsigned long int n = 0; n < totSize; n++)
            {
                locData[n] /= (T)number;
            }
            return *this;
        }

        template<typename T2>
        Tensor<T, Rank>& operator*=(const T2 number)
        {
            if (allocated)
            for(unsigned long int n = 0; n < totSize; n++)
            {
                locData[n] *= (T)number;
            }
            return *this;
        }

        Tensor<T, Rank> operator+(const Tensor<T, Rank>& locTensor) const
        {
            if (locTensor.Dimensions == Dimensions)
            {
                Tensor<T,Rank> myReturn(locTensor);

                for(unsigned long int n = 0; n < totSize; n++)
                {
                    myReturn[n] += locData[n];
                }
                return myReturn;
            }
            else
            {
                std::cout << "ERROR: Tensor<T," << Rank << ">: Different tensor size in summation operator!"
                        << std::endl;
                exit(133);
                return *this;
            }
        }
        Tensor<T, Rank> operator-(const Tensor<T, Rank>& locTensor) const
        {
            if (locTensor.Dimensions == Dimensions)
            {
                Tensor<T,Rank> myReturn(locTensor);

                for(unsigned long int n = 0; n < totSize; n++)
                {
                    myReturn[n] -= locData[n];
                }
                return myReturn;
            }
            else
            {
                std::cout << "ERROR: Tensor<T," << Rank << ">: Different tensor size in subtraction operator!"
                        << std::endl;
                exit(133);
                return *this;
            }
        }
        Tensor<T, Rank> operator*(double val) const
        {
            Tensor<T, Rank> myReturn(*this);
            for(unsigned long int n = 0; n < totSize; n++)
            {
                myReturn[n] *= val;
            }
            return myReturn;
        }
        Tensor<T, Rank> operator/(double val) const
        {
            Tensor<T, Rank> myReturn(*this);
            for(unsigned long int n = 0; n < totSize; n++)
            {
                myReturn[n] /= val;
            }
            return myReturn;
        }
        bool IsAllocated() const
        {
            return allocated;
        }
    protected:
        std::array<long int, Rank> Dimensions;
        unsigned long int totSize;
        bool allocated;
        T*  locData;
    private:
        long int Index(const std::initializer_list<long int> position) const
        {
            long int locIndex = *(position.begin());
            for(long int n = 1; n < Rank; n++)
            {
                locIndex *= Dimensions[n];
                locIndex += *(position.begin() + n);
            }
            return locIndex;
        }
    };

}
#endif