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

#ifndef NODEA_H
#define NODEA_H

#include "Tools/Includes.h"
namespace opensim
{
/***************************************************************/

struct AdvEntry                                                                 /// Structure for storing the individual phase field derivatives entry. Used in the NodeA class as a storage unit.
{
    int    index;                                                               ///  Phase field index
    double dXY;                                                                 ///  Second derivative w.r.t. X and Y
    double dXZ;                                                                 ///  Second derivative w.r.t. X and Z
    double dYZ;                                                                 ///  Second derivative w.r.t. Y and Z
    double dXYZ;                                                                ///  Third derivative w.r.t. X, Y and z
};

/***************************************************************/
class NodeA                                                                     /// Stores the phase fields derivatives at a grid point. Provide access and manipulation methods for the entries.
{
 public:
    double operator[](const int    n);                                          /// Index operator for accessing the n's phase field derivatives

    void   set_dXY (const int n, const double value);                           /// Set second derivative w.r.t X and Y
    void   set_dXZ (const int n, const double value);                           /// Set second derivative w.r.t X and Z
    void   set_dYZ (const int n, const double value);                           /// Set second derivative w.r.t Y and Z
    void   set_dXYZ(const int n, const double value);                           /// Set third derivative w.r.t X, Y and Z

    void   set(const int n, const double dXY,
                            const double dXZ,
                            const double dYZ,
                            const double dXYZ);                                 /// Set all derivatives simultaneously.

    void   add_dXY (const int n, const double value);                           /// Increment second derivative w.r.t X and Y
    void   add_dXZ (const int n, const double value);                           /// Increment second derivative w.r.t X and Z
    void   add_dYZ (const int n, const double value);                           /// Increment second derivative w.r.t Y and Z
    void   add_dXYZ(const int n, const double value);                           /// Increment third derivative w.r.t X, Y and Z

    double get_dXY (const int n);                                               /// Return second derivative w.r.t X and Y
    double get_dXZ (const int n);                                               /// Return second derivative w.r.t X and Z
    double get_dYZ (const int n);                                               /// Return second derivative w.r.t Y and Z
    double get_dXYZ(const int n);                                               /// Return second derivative w.r.t X, Y and Z

    void   clear()                                                              /// Emptys the storage.
           {Fields.clear();};
    int    size()                                                               /// Returns the size of storage.
           {return Fields.size();};
    typedef std::vector<AdvEntry>::iterator iterator;                           /// Iterator over phase field derivatives storage
    typedef std::vector<AdvEntry>::const_iterator citerator;                    /// Constant iterator over phase field derivatives storage
    iterator begin()
             {return Fields.begin();};                                          /// Iterator to the begin of phase field derivatives storage
    iterator end()                                                              /// Iterator to the end of phase field derivatives storage
             {return Fields.end();};
    citerator cbegin()
             {return Fields.begin();};                                          /// Constant iterator to the begin of phase field derivatives storage
    citerator cend()                                                            /// Constant iterator to the end of phase field derivatives storage
             {return Fields.end();};

 protected:
 private:
    std::vector<AdvEntry> Fields;                                               /// List of nonvanishing phase field derivatives.
};

inline void NodeA::set_dXY(const int n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dXY = value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index = n;
    NewEntry.dXY   = value;
    NewEntry.dXZ   = 0.0;
    NewEntry.dYZ   = 0.0;
    NewEntry.dXYZ  = 0.0;
    Fields.push_back(NewEntry);
}

inline void NodeA::set_dXZ(const int n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dXZ = value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index = n;
    NewEntry.dXZ   = value;
    NewEntry.dXY   = 0.0;
    NewEntry.dYZ   = 0.0;
    NewEntry.dXYZ  = 0.0;
    Fields.push_back(NewEntry);
}

inline void NodeA::set_dYZ(const int n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dYZ = value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index = n;
    NewEntry.dYZ   = value;
    NewEntry.dXY   = 0.0;
    NewEntry.dXZ   = 0.0;
    NewEntry.dXYZ  = 0.0;
    Fields.push_back(NewEntry);
}

inline void NodeA::set_dXYZ(const int n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dXYZ = value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index = n;
    NewEntry.dXYZ  = value;
    NewEntry.dXY   = 0.0;
    NewEntry.dXZ   = 0.0;
    NewEntry.dYZ   = 0.0;
    Fields.push_back(NewEntry);
}

inline void NodeA::set(const int n, const double dXY,
                             const double dXZ,
                             const double dYZ,
                             const double dXYZ)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dXY  = dXY;
        i->dXZ  = dXZ;
        i->dYZ  = dYZ;
        i->dXYZ = dXYZ;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index     = n;
    NewEntry.dXY   = dXY;
    NewEntry.dXZ   = dXZ;
    NewEntry.dYZ   = dYZ;
    NewEntry.dXYZ  = dXYZ;
    Fields.push_back(NewEntry);
}

inline void NodeA::add_dXY(const int n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dXY += value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index     = n;
    NewEntry.dXY   = value;
    NewEntry.dXZ   = 0.0;
    NewEntry.dYZ   = 0.0;
    NewEntry.dXYZ  = 0.0;
    Fields.push_back(NewEntry);
}

inline void NodeA::add_dXZ(const int n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dXZ += value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index = n;
    NewEntry.dXZ   = value;
    NewEntry.dXY   = 0.0;
    NewEntry.dYZ   = 0.0;
    NewEntry.dXYZ  = 0.0;
    Fields.push_back(NewEntry);
}

inline void NodeA::add_dYZ(const int n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dYZ += value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index = n;
    NewEntry.dYZ   = value;
    NewEntry.dXY   = 0.0;
    NewEntry.dXZ   = 0.0;
    NewEntry.dXYZ  = 0.0;
    Fields.push_back(NewEntry);
}

inline void NodeA::add_dXYZ(const int n, const double value)
{
    for (iterator i = Fields.begin(); i < Fields.end(); ++i)
    if (i->index == n)
    {
        i->dXYZ += value;
        return;
    }

    AdvEntry NewEntry;
    NewEntry.index = n;
    NewEntry.dXYZ  = value;
    NewEntry.dXY   = 0.0;
    NewEntry.dXZ   = 0.0;
    NewEntry.dYZ   = 0.0;
    Fields.push_back(NewEntry);
}

inline double NodeA::get_dXY(const int n)
{
    for (citerator i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if (i->index == n) return i->dXY;
    }
    return 0.0;
}

inline double NodeA::get_dXZ(const int n)
{
    for (citerator i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if (i->index == n) return i->dXZ;
    }
    return 0.0;
}

inline double NodeA::get_dYZ(const int n)
{
    for (citerator i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if (i->index == n) return i->dYZ;
    }
    return 0.0;
}

inline double NodeA::get_dXYZ(const int n)
{
    for (citerator i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if (i->index == n) return i->dXYZ;
    }
    return 0;
}
}// namespace openphase
#endif
