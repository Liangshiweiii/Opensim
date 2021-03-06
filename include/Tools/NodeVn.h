

#ifndef NODEVN_H
#define NODEVN_H

#include "Tools/Includes.h"
namespace opensim
{

template<int N>
struct VectorVnEntry                                                            /// Individual entries. Used in the NodeVn class as a storage unit.
{
    double value[N];
    int index;
};

/***************************************************************/
template<int N>
class NodeVn                                                                    /// Basic class to store the vector valued quantities for each phase field. Provide access and manipulation methods for the entrys.
{
 public:
    void   set      (const int idx, const int ss, const double value);          /// Set component ss for phase field index idx to value.
    void   set      (const int idx, const dVector<N> value);                    /// Set dVector<N> for phase field index idx.
    void   set      (const int idx, const double value);                        /// Set all entries to single value for phase field index idx.
    void   set_to_zero (const int idx);                                         /// Set all entries to zero for phase field index idx.
    double get      (const int idx,  const int ss) const;                       /// Return value of component ss for phase field index idx.
    dVector<N> get      (const int idx) const;                                  /// Return value of component ss for phase field index idx.
    void   add    (const int n, const double value[N]);                         /// Increment value using two indeces.
    void   change_index  (const int n, const int m);

    NodeVn&   operator=(const NodeVn& n);
    NodeVn   operator+(const  NodeVn& n) const;                                 /// Plus operator. Takes as input another Node type entry.
    NodeVn   operator-(const  NodeVn& n) const;                                 /// Minus operator. Takes as input another Node type entry.
    NodeVn   operator*(const double n) const;                                   /// Multiply all fields by a number.

    NodeVn   add     (const NodeVn& value) const;                               /// Add two nodes.

    void   clear(){VectorVnFields.clear();};                                    /// Empties the vector fields.
    int    size() const {return VectorVnFields.size();};                        /// Returns the size of vector fields.
    typedef typename std::vector<VectorVnEntry<N>>::iterator iterator;          /// Iterator over the vector fields
    typedef typename std::vector<VectorVnEntry<N>>::const_iterator citerator;   /// Constant iterator over the vector fields
    iterator begin() {return VectorVnFields.begin();};                          /// Iterator to the begin of vector fields
    iterator end()   {return VectorVnFields.end();};                            /// Iterator to the end of vector fields
    citerator cbegin() const {return VectorVnFields.begin();};                  /// Constant iterator to the begin of vector fields
    citerator cend()   const {return VectorVnFields.end();};                    /// Constant iterator to the end of vector fields

 protected:
 private:
    std::vector<VectorVnEntry<N>> VectorVnFields;                               /// List of nonvanishing vector fields.
};

template<int N>
inline void NodeVn<N>::set(const int idx, const int ss, const double value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == idx)
    {
        i->value[ss] = value;
        return;
    }

    VectorVnEntry<N> NewEntry;
    NewEntry.index = idx;
    for(int ii = 0; ii < N; ii++) NewEntry.value[ii] = 0.0;
    NewEntry.value[ss] = value;

    VectorVnFields.push_back(NewEntry);
}

template<int N>
inline void NodeVn<N>::set(const int idx, const dVector<N> value)
{
    for (iterator i = begin(); i < end(); ++i)
    if (i->index == idx)
    {
        for (unsigned int ss = 0; ss < N; ss++)
        {
            i->value[ss] = value[ss];
        }
        return;
    }

    VectorVnEntry<N> NewEntry;
    NewEntry.index = idx;
    for(int ii = 0; ii < N; ii++)
    {
        NewEntry.value[ii] = value[ii];
    }
    VectorVnFields.push_back(NewEntry);
}

template<int N>
inline void NodeVn<N>::set(const int idx, const double value)
{
    for (iterator i = begin(); i < end(); ++i)
    {
        if (i->index == idx)
        {
            for (int ii = 0; ii < N; ++ii)
            {
                i->value[ii] = value;
            }
            return;
        }
    }

    VectorVnEntry<N> NewEntry;
    NewEntry.index = idx;
    for(int ii = 0; ii < N; ii++) NewEntry.value[ii] = value;

    VectorVnFields.push_back(NewEntry);
}

template<int N>
inline void NodeVn<N>::set_to_zero(const int idx)
{
    for (iterator i = begin(); i < end(); ++i)
    {
        if (i->index == idx)
        {
            for (int ii = 0; ii < N; ++ii)
            {
                i->value[ii] = 0.0;
            }
            return;
        }
    }

    VectorVnEntry<N> NewEntry;
    NewEntry.index = idx;
    for(int ii = 0; ii < N; ii++) NewEntry.value[ii] = 0.0;

    VectorVnFields.push_back(NewEntry);
}

template<int N>
inline void NodeVn<N>::change_index(const int n, const int m)
{
    // Set values and index of m to index n
    // Warning: Values are averaged arithmetically!
    if (n != m)
    {
        double returndV[N];
        for (unsigned int ss = 0; ss < N; ss++) returndV[ss] = 0.0;

        bool found = false;
        for (auto i = cbegin(); i < cend(); ++i)
        {
            if (i->index == m)
            {
                for (unsigned int ss = 0; ss < N; ss++)
                {
                    returndV[ss] = (i->value)[ss];
                }
                found = true;
                break;
            }
        }
        if (found == true)
        {
            for (auto i = begin(); i < end(); ++i)
            if (i->index == m)
            {
                for (unsigned int ss = 0; ss < N; ss++)
                {
                    (i->value)[ss] = 0.5*((i->value)[ss]+returndV[ss]);
                }
                i->index = n;
                break;
            }
        }
    }
}

template<int N>
inline double NodeVn<N>::get(const int idx, const int ss) const
{
    for (citerator i = cbegin(); i < cend(); ++i)
    if (i->index == idx)
    {
        return i->value[ss];
    }
    return 0.0;
}

template<int N>
inline dVector<N> NodeVn<N>::get(const int idx) const
{
    dVector<N> returndV; returndV.set_to_zero();
    for (citerator i = cbegin(); i < cend(); ++i)
    {
        if (i->index == idx)
        {
            for (unsigned int ss = 0; ss < N; ss++)
            {
                returndV[ss] = (i->value)[ss];
            }
            break;
        }
    }
    return returndV;
}

template<int N>
inline void NodeVn<N>::add(const int n, const double summand[N])
{
    for (auto i = begin(); i < end(); ++i)
    if (i->index == n)
    {
        for(int ss = 0; ss < N; ++ss)
        {
            i->value[ss] += summand[ss];
        }
        return;
    }

    VectorVnEntry<N> NewEntry;
    NewEntry.index = n;
    for (int ss = 0; ss < N; ss ++) NewEntry.value[ss] = summand[ss];

    VectorVnFields.push_back(NewEntry);
}

template<int N>
inline NodeVn<N> NodeVn<N>::add(const NodeVn& n) const
{
    NodeVn<N> result = n;

    for (auto i = cbegin(); i < cend(); ++i)
    {
        result.add(i->index, i->value);
    }
    return result;
}

template<int N>
inline NodeVn<N>& NodeVn<N>::operator=(const NodeVn& n)
{
    VectorVnFields = n.VectorVnFields;
    return *this;
}

template<int N>
inline NodeVn<N> NodeVn<N>::operator+(const NodeVn& n) const
{
    NodeVn<N> result = n;

    for (auto i = cbegin(); i < cend(); ++i)
    {
        result.add(i->index, i->value);
    }
    return result;
}

template<int N>
inline NodeVn<N> NodeVn<N>::operator-(const NodeVn& n) const
{
    NodeVn<N> result = n;

    for (auto i = result.begin(); i < result.end(); ++i)
    for(int ss = 0; ss < N; ++ss)
    {
        i->value[ss] = -i->value[ss];
    }

    for (auto i = cbegin(); i < cend(); ++i)
    {
        result.add(i->index, i->value);
    }
    return result;
}

template<int N>
inline NodeVn<N> NodeVn<N>::operator*(const double n) const
{
    NodeVn<N> result = *this;

    for(auto i = result.begin(); i < result.end(); ++i)
    for(int ss = 0; ss < N; ++ss)
    {
        i->value[ss] *= n;
    }

    return result;
}

} //namespace opensim
#endif
