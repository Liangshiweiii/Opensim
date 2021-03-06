/* This file is part of OPENPHASE BASE.
 * Copyright (c) 2020 OpenPhase Solutions GmbH, Bochum, Germany
 * For more details visit https://www.openphase-solutions.com
 * 
 * OPENPHASE BASE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *  
 * OPENPHASE BASE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *  
 * You should have received a copy of the GNU General Public License
 * along with OPENPHASE BASE.  If not, see <http://www.gnu.org/licenses/>. */

#ifndef STORAGE_H
#define STORAGE_H

namespace opensim
{

template <class T>
class Storage  /// 1D storage template class. Can handle any type of values
{
 public:
    Storage()
    {
        Array = nullptr;
    }
    T& operator[](const int x)
    {
#ifdef DEBUG
        if(x > Size_X)
        {
            std::cout << "Error in Storage::operator[]\n"
                 << "Access beyond storage range. x = " << x
                 << " > Size_X = " << Size_X <<"\nTerminating!!!"<<std::endl;
            exit(13);
        }
#endif
        return Array[x];
    }
    T const& operator[](const int x) const
    {
#ifdef DEBUG
        if(x > Size_X)
        {
            std::cout << "Error in Storage::operator[]\n"
                 << "Access beyond storage range. x = " << x
                 << " > Size_X = " << Size_X <<"\nTerminating!!!"<<std::endl;
            exit(13);
        }
#endif
        return Array[x];
    }
    void Allocate(const int nx)
    {
        Size_X = nx;
        Array = new T[Size_X] ();
    }
    void Reallocate(const int nx)
    {
        delete[] Array;
        Size_X = nx;
        Array = new T[Size_X] ();
    }
    int Size() const
    {
        return Size_X;
    }
    bool IsNotAllocated() const
    {
        return (Array == nullptr);
    }
    bool IsAllocated() const
    {
        return !(Array == nullptr);
    }
    ~Storage()
    {
        delete[] Array;
    }

 protected:
 private:
    T* Array;
    int Size_X;
};

}// namespace opensim
#endif
