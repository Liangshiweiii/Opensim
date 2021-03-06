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


#ifndef NODEV_H
#define NODEV_H

#include "Tools/Includes.h"
namespace opensim
{

struct VectorEntry                                                              /// Individual vector entry. Used in the NodeV class as a storage unit.
{
    double X;                                                                   /// X component of the vector
    double Y;                                                                   /// Y component of the vector
    double Z;                                                                   /// Z component of the vector
    union                                                                       /// First index
    {
        int index;
        int indexA;
    };
    int     indexB;                                                             /// Second index
};

/***************************************************************/

class NodeV                                                                     /// Basic class to store the vector valued quantities for each phase field. Provide access and manipulation methods for the entrys.
{
 public:
    NodeV  operator+(const  NodeV& n) const;                                    /// Plus operator. Takes as input another NodeV type entry. Sums only x, y and z components.
    NodeV  operator-(const  NodeV& n) const;                                    /// Minus operator. Takes as input another NodeV type entry. Subtructs only x, y and z components.

    NodeV  Xreflected() const;                                                  /// Returns the NodeV with the sign of X components of all vectors changed
    NodeV  Yreflected() const;                                                  /// Returns the NodeV with the sign of Y components of all vectors changed
    NodeV  Zreflected() const;                                                  /// Returns the NodeV with the sign of Z components of all vectors changed
    //NodeV  operator*(const  double n) const;                                  /// Multiply operator. Takes as input the multiplier. Multiplies only x, y and z components
    //NodeV  operator/(const  double n) const;                                  /// Divide operator. Takes as input the divizor. Divides only x, y and z components

    void   set    (const int n, const dVector3);                                /// Set all components using dVector3 type
    void   set_X    (const int n, const double value);                          /// Set X component.
    void   set_Y    (const int n, const double value);                          /// Set Y component.
    void   set_Z    (const int n, const double value);                          /// Set Z component.
    void   set      (const int n, const double valueX,
                                  const double valueY,
                                  const double valueZ);                         /// Set all components.

    dVector3 get    (const int n) const;                                        /// Returns all components as dVector3 type
    double get_X    (const int n) const;                                        /// Return X component.
    double get_Y    (const int n) const;                                        /// Return Y component.
    double get_Z    (const int n) const;                                        /// Return Z component.
    void   get      (const int n, double& valueX,
                                  double& valueY,
                                  double& valueZ) const;                        /// Return all components.

    void   add_X    (const int n, const double value);                          /// Add to X component.
    void   add_Y    (const int n, const double value);                          /// Add to Y component.
    void   add_Z    (const int n, const double value);                          /// Add to Z component.
    void   add      (const int n, const double valueX,
                                  const double valueY,
                                  const double valueZ);                         /// Add to all components.

    /// For dual indeces we assume antisymmetric properties f(n,m) = -f(m,n).

    void   set    (const int n, const int m, const dVector3);                   /// Set all components using dVector3 type
    void   set_X    (const int n, const int m, const double value);             /// Set X component.
    void   set_Y    (const int n, const int m, const double value);             /// Set Y component.
    void   set_Z    (const int n, const int m, const double value);             /// Set Z component.
    void   set      (const int n, const int m, const double valueX,
                                               const double valueY,
                                               const double valueZ);            /// Set all components.

    dVector3 get    (const int n, const int m) const;                           /// Returns all components as dVector3 type
    double get_X    (const int n, const int m) const;                           /// Return X component.
    double get_Y    (const int n, const int m) const;                           /// Return Y component.
    double get_Z    (const int n, const int m) const;                           /// Return Z component.
    void   get      (const int n, const int m, double& valueX,
                                               double& valueY,
                                               double& valueZ) const;           /// Return all components.

    void   add_X    (const int n, const int m, const double value);             /// Add to X component.
    void   add_Y    (const int n, const int m, const double value);             /// Add to Y component.
    void   add_Z    (const int n, const int m, const double value);             /// Add to Z component.
    void   add      (const int n, const int m, const double valueX,
                                               const double valueY,
                                               const double valueZ);            /// Add to all components.

    void   clear(){VectorFields.clear();};                                      /// Empties the vector fields.
    int    size() const {return VectorFields.size();};                          /// Returns the size of vector fields.
    typedef std::vector<VectorEntry>::iterator iterator;                        /// Iterator over the vector fields
    typedef std::vector<VectorEntry>::const_iterator citerator;                 /// Constant iterator over the vector fields
    iterator begin() {return VectorFields.begin();};                            /// Iterator to the begin of vector fields
    iterator end()   {return VectorFields.end();};                              /// Iterator to the end of vector fields
    citerator cbegin() const {return VectorFields.cbegin();};                   /// Constant iterator to the begin of vector fields
    citerator cend() const   {return VectorFields.cend();};                     /// Constant iterator to the end of vector fields

 protected:
 private:
    std::vector<VectorEntry> VectorFields;                                      /// List of nonvanishing vector fields.
};

inline NodeV NodeV::operator+(const  NodeV& n) const
{
    NodeV result = n;
    for (citerator i = VectorFields.cbegin(); i < VectorFields.cend(); ++i)
    {
        result.add(i->indexA, i->indexB, i->X, i->Y, i->Z);
    }
    return result;
}

inline NodeV NodeV::operator-(const NodeV& n) const
{
    NodeV result = n;

    for (iterator i = result.begin(); i < result.end(); ++i)
    {
        i->X = -i->X;
        i->Y = -i->Y;
        i->Z = -i->Z;
    }

    for (citerator i = VectorFields.cbegin(); i < VectorFields.cend(); ++i)
    {
        result.add(i->indexA, i->indexB, i->X, i->Y, i->Z);
    }
    return result;
}

inline NodeV NodeV::Xreflected() const
{
    NodeV result = *this;

    for (iterator i = result.begin(); i < result.end(); ++i)
    {
        i->X = -i->X;
    }
    return result;
}

inline NodeV NodeV::Yreflected() const
{
    NodeV result = *this;

    for (iterator i = result.begin(); i < result.end(); ++i)
    {
        i->Y = -i->Y;
    }
    return result;
}

inline NodeV NodeV::Zreflected() const
{
    NodeV result = *this;

    for (iterator i = result.begin(); i < result.end(); ++i)
    {
        i->Z = -i->Z;
    }

    return result;
}

inline void NodeV::set(const int n, const dVector3 value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == n)
    {
        i->X = value[0];
        i->Y = value[1];
        i->Z = value[2];
        return;
    }

    VectorEntry NewEntry;
    NewEntry.index  = n;
    NewEntry.indexB = 0;
    NewEntry.X     = value[0];
    NewEntry.Y     = value[1];
    NewEntry.Z     = value[2];

    VectorFields.push_back(NewEntry);
}

inline void NodeV::set_X(const int n, const double value)
{
    for (iterator i =  begin(); i < end(); ++i)
    if (i->index == n)
    {
        i->X = value;
        return;
    }

    VectorEntry NewEntry;
    NewEntry.index = n;
    NewEntry.indexB = 0;
    NewEntry.X     = value;
    NewEntry.Y     = 0.0;
    NewEntry.Z     = 0.0;
    VectorFields.push_back(NewEntry);
}

inline void NodeV::set_Y(const int n, const double value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == n)
    {
        i->Y = value;
        return;
    }

    VectorEntry NewEntry;
    NewEntry.index = n;
    NewEntry.indexB = 0;
    NewEntry.X     = 0.0;
    NewEntry.Y     = value;
    NewEntry.Z     = 0.0;
    VectorFields.push_back(NewEntry);
}

inline void NodeV::set_Z(const int n, const double value)
{
    for (iterator i =  begin(); i < end(); ++i)
    if (i->index == n)
    {
        i->Z = value;
        return;
    }

    VectorEntry NewEntry;
    NewEntry.index = n;
    NewEntry.indexB = 0;
    NewEntry.X     = 0.0;
    NewEntry.Y     = 0.0;
    NewEntry.Z     = value;
    VectorFields.push_back(NewEntry);
}

inline void NodeV::set(const int n, const double valueX,
                                    const double valueY,
                                    const double valueZ)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == n)
    {
        i->X = valueX;
        i->Y = valueY;
        i->Z = valueZ;
        return;
    }

    VectorEntry NewEntry;
    NewEntry.index = n;
    NewEntry.indexB = 0;
    NewEntry.X     = valueX;
    NewEntry.Y     = valueY;
    NewEntry.Z     = valueZ;

    VectorFields.push_back(NewEntry);
}

inline dVector3 NodeV::get(const int n) const
{
    dVector3 result;
    result.set_to_zero();

    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->index == n)
    {

        result[0] = i->X;
        result[1] = i->Y;
        result[2] = i->Z;
    }
    return result;
}

inline double NodeV::get_X(const int n) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->index == n) return i->X;

    return 0.0;
}

inline double NodeV::get_Y(const int n) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->index == n) return i->Y;

    return 0.0;
}

inline double NodeV::get_Z(const int n) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->index == n) return i->Z;

    return 0.0;
}

inline void NodeV::get(const int n, double& valueX, double& valueY,
        double& valueZ) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->index == n)
    {
        valueX = i->X;
        valueY = i->Y;
        valueZ = i->Z;
        return;
    }

    valueX = 0.0;
    valueY = 0.0;
    valueZ = 0.0;
}

inline void NodeV::add_X(const int n, const double value)
{
    for (iterator i =  begin(); i < end(); ++i)
    if (i->index == n)
    {
        i->X += value;
        return;
    }

    VectorEntry NewEntry;
    NewEntry.index = n;
    NewEntry.indexB = 0;
    NewEntry.X     = value;
    NewEntry.Y     = 0.0;
    NewEntry.Z     = 0.0;
    VectorFields.push_back(NewEntry);
}

inline void NodeV::add_Y(const int n, const double value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == n)
    {
        i->Y += value;
        return;
    }

    VectorEntry NewEntry;
    NewEntry.index = n;
    NewEntry.indexB = 0;
    NewEntry.X     = 0.0;
    NewEntry.Y     = value;
    NewEntry.Z     = 0.0;
    VectorFields.push_back(NewEntry);
}

inline void NodeV::add_Z(const int n, const double value)
{
    for (iterator i =  begin(); i < end(); ++i)
    if (i->index == n)
    {
        i->Z += value;
        return;
    }

    VectorEntry NewEntry;
    NewEntry.index = n;
    NewEntry.indexB = 0;
    NewEntry.X     = 0.0;
    NewEntry.Y     = 0.0;
    NewEntry.Z     = value;
    VectorFields.push_back(NewEntry);
}

inline void NodeV::add(const int n, const double valueX,
                                    const double valueY,
                                    const double valueZ)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == n)
    {
        i->X     += valueX;
        i->Y     += valueY;
        i->Z     += valueZ;
        return;
    }

    VectorEntry NewEntry;
    NewEntry.index = n;
    NewEntry.indexB = 0;
    NewEntry.X     = valueX;
    NewEntry.Y     = valueY;
    NewEntry.Z     = valueZ;

    VectorFields.push_back(NewEntry);
}

inline void NodeV::set(const int n, const int m, const dVector3 value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->X = value[0];
        i->Y = value[1];
        i->Z = value[2];
        return;
    }
    else if (i->indexA == m && i->indexB == n)
    {
        i->X = -value[0];
        i->Y = -value[1];
        i->Z = -value[2];
        return;
    }

    VectorEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.X      = value[0];
    NewEntry.Y      = value[1];
    NewEntry.Z      = value[2];

    VectorFields.push_back(NewEntry);
}

inline void NodeV::set_X(const int n, const int m, const double value)
{
    for (iterator i =  begin(); i < end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->X = value;
        return;
    }
    else if (i->indexA == m && i->indexB == n)
    {
        i->X = -value;
        return;
    }
    VectorEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.X     = value;
    NewEntry.Y     = 0.0;
    NewEntry.Z     = 0.0;
    VectorFields.push_back(NewEntry);
}

inline void NodeV::set_Y(const int n, const int m, const double value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->Y = value;
        return;
    }
    else if (i->indexA == m && i->indexB == n)
    {
        i->Y = -value;
        return;
    }

    VectorEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.X     = 0.0;
    NewEntry.Y     = value;
    NewEntry.Z     = 0.0;
    VectorFields.push_back(NewEntry);
}

inline void NodeV::set_Z(const int n, const int m, const double value)
{
    for (iterator i =  begin(); i < end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->Z = value;
        return;
    }
    else if (i->indexA == m && i->indexB == n)
    {
        i->Z = -value;
        return;
    }

    VectorEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.X     = 0.0;
    NewEntry.Y     = 0.0;
    NewEntry.Z     = value;
    VectorFields.push_back(NewEntry);
}

inline void NodeV::set(const int n, const int m, const double valueX,
                                                 const double valueY,
                                                 const double valueZ)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->X = valueX;
        i->Y = valueY;
        i->Z = valueZ;
        return;
    }
    else if (i->indexA == m && i->indexB == n)
    {
        i->X = -valueX;
        i->Y = -valueY;
        i->Z = -valueZ;
        return;
    }

    VectorEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.X     = valueX;
    NewEntry.Y     = valueY;
    NewEntry.Z     = valueZ;

    VectorFields.push_back(NewEntry);
}

inline dVector3 NodeV::get(const int n, const int m) const
{
    dVector3 result;
    result.set_to_zero();

    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->indexA == n && i->indexB == m)
    {

        result[0] = i->X;
        result[1] = i->Y;
        result[2] = i->Z;
    }
    else if (i->indexA == m && i->indexB == n)
    {
        dVector3 result;
        result[0] = -i->X;
        result[1] = -i->Y;
        result[2] = -i->Z;
    }

    return result;
}

inline double NodeV::get_X(const int n, const int m) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->indexA == n && i->indexB == m) return i->X;
    else if (i->indexA == m && i->indexB == n) return -i->X;

    return 0.0;
}

inline double NodeV::get_Y(const int n, const int m) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->indexA == n && i->indexB == m) return i->Y;
    else if (i->indexA == m && i->indexB == n) return -i->Y;

    return 0.0;
}

inline double NodeV::get_Z(const int n, const int m) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->indexA == n && i->indexB == m) return i->Z;
    else if (i->indexA == m && i->indexB == n) return -i->Z;

    return 0.0;
}

inline void NodeV::get(const int n, const int m, double& valueX, double& valueY,
        double& valueZ) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        valueX = i->X;
        valueY = i->Y;
        valueZ = i->Z;
        return;
    }
    else if (i->indexA == m && i->indexB == n)
    {
        valueX = -i->X;
        valueY = -i->Y;
        valueZ = -i->Z;
        return;
    }
    valueX = 0.0;
    valueY = 0.0;
    valueZ = 0.0;
}

inline void NodeV::add_X(const int n, const int m, const double value)
{
    for (iterator i =  begin(); i < end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->X += value;
        return;
    }
    else if (i->indexA == m && i->indexB == n)
    {
        i->X -= value;
        return;
    }

    VectorEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.X      = value;
    NewEntry.Y      = 0.0;
    NewEntry.Z      = 0.0;
    VectorFields.push_back(NewEntry);
}

inline void NodeV::add_Y(const int n, const int m, const double value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->Y += value;
        return;
    }
    else if (i->indexA == m && i->indexB == n)
    {
        i->Y -= value;
        return;
    }

    VectorEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.X      = 0.0;
    NewEntry.Y      = value;
    NewEntry.Z      = 0.0;
    VectorFields.push_back(NewEntry);
}

inline void NodeV::add_Z(const int n, const int m, const double value)
{
    for (iterator i =  begin(); i < end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->Z += value;
        return;
    }
    else if (i->indexA == m && i->indexB == n)
    {
        i->Z -= value;
        return;
    }

    VectorEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.X      = 0.0;
    NewEntry.Y      = 0.0;
    NewEntry.Z      = value;
    VectorFields.push_back(NewEntry);
}

inline void NodeV::add(const int n, const int m, const double valueX,
                                                 const double valueY,
                                                 const double valueZ)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->X += valueX;
        i->Y += valueY;
        i->Z += valueZ;
        return;
    }
    else if (i->indexA == m && i->indexB == n)
    {
        i->X -= valueX;
        i->Y -= valueY;
        i->Z -= valueZ;
        return;
    }
    VectorEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.X      = valueX;
    NewEntry.Y      = valueY;
    NewEntry.Z      = valueZ;

    VectorFields.push_back(NewEntry);
}

} //namespace openphase
#endif
