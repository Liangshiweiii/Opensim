#ifndef NODE_H
#define NODE_H

#include "Tools/Includes.h"

namespace opensim
{
/**********************************************************/

struct FieldEntry                                                               ///< Structure for storing the field entry. Used in the Node class as a storage unit.
{
    union                                                                       ///< First index
    {
        int index;
        int indexA;
    };
    union                                                                       ///< Second index
    {
        int indexP;
        int indexB;
    };
    union                                                                       ///< First value to store
    {
        double value;
        double value1;
    };
    union                                                                       ///< Second value to store
    {
        double laplacian;
        double value2;
    };
};

/**********************************************************/
class Node                                                                      ///< 在网格点中储存相场的数值以及它们的倒数等等，提供了一个很好的途径和方法来处理相场的输入
{
 public:
    Node()                                                                      ///<构造函数，用于为至少3个字段分配空间
    {
        Fields.reserve(3);
        flag = 0;
    }
    ~Node()
    {

    }
    double operator[](const int n) const;                                       ///< 用于访问n的相场值的索引操作符

    Node   operator+(const  Node& n) const;                                     ///< Plus operator. Takes as input another Node type entry.
    Node   operator-(const  Node& n) const;                                     ///< Minus operator. Takes as input another Node type entry.
    Node   operator*(const double n) const;                                     ///< Multiply all fields by a number.
    Node&  operator=(const Node& n);                                            ///< Assignement operator
    int    finalize(void);                                                      ///< Used for phase fields. Adjusts phase field values to (0, 1] interval.

    void   set    (const int n, const int m, const double value);               ///< Set value using two indeces.
    double get    (const int n, const int m) const;                             ///< Return value using two indeces.
    void   add    (const int n, const int m, const double value);               ///< Increment value using two indeces.

    void   set2    (const int n, const int m, const double value2);             ///< Set value2 using two indeces.
    double get2    (const int n, const int m) const;                            ///< Return value2 using two indeces.
    void   add2    (const int n, const int m, const double value2);             ///< Increment value2 using two indeces.

    void   set    (const int n, const double value);                            ///< Set value using single index.
    double get    (const int n) const;                                          ///< Return value using single index.
    void   add    (const int n, const double value);                            ///< Increment value using single index.

    void   set2    (const int n, const double value2);                          ///< Set value2 using single index.
    double get2    (const int n) const;                                         ///< Return value2 using single index.
    void   add2    (const int n, const double value2);                          ///< Increment value2 using single index.

    void   set_sym(const int n, const int m, const double value);               ///< Set value in symmetric case: f(n,m) = f(m,n).
    void   set_sym2(const int n, const int m, const double value);               ///< Set value2 in symmetric case: f(n,m) = f(m,n).
    double get_sym(const int n, const int m) const;                             ///< Return value in symmetric case: f(n,m) = f(m,n).
    double get_sym2(const int n, const int m) const;                            ///< Return value2 in symmetric case: f(n,m) = f(m,n).
    void   add_sym(const int n, const int m,  const double value);              ///< Increment value in symmetric case: f(n,m) = f(m,n).

    void   set_asym(const int n, const int m, const double value);              ///< Set value in antisymmetric case: f(n,m) = -f(m,n).
    double get_asym(const int n, const int m) const;                            ///< Return value in antisymmetric case: f(n,m) = -f(m,n).
    double get_asym2(const int n, const int m) const;                            ///< Return value2 in antisymmetric case: f(n,m) = -f(m,n).
    void   add_asym(const int n, const int m,  const double value);             ///< Increment value in antisymmetric case: f(n,m) = -f(m,n).

    FieldEntry get_max(void) const;                                             ///< Returns vector of FieldEntry with max positive abs value (value1)

    Node   add     (const Node& value) const;                                   ///< Add two nodes.
    Node   add_sym (const Node& value) const;                                   ///< Add two nodes in symmetric case: f(n,m) = f(m,n).
    Node   add_asym(const Node& value) const;                                   ///< Add two nodes in antisymmetric case: f(n,m) = -f(m,n).
    Node   add_exist(const Node& value) const;                                  ///< Add only entries existing in two nodes simultaneously.
    Node   add_sym_exist(const Node& value) const;                              ///< Add only entries existing in two nodes simultaneously in symmetric case: f(n,m) = f(m,n).
    Node   add_asym_exist(const Node& value) const;                             ///< Add only entries existing in two nodes simultaneously in antisymmetric case: f(n,m) = -f(m,n).
    void   clear()                                                              ///< Emptys the field storage. Sets flag to 0.
    {
        Fields.clear();
        flag = 0;
    };

    unsigned int size() const                                                   ///< Returns the size of storage.
    {
        return Fields.size();
    };
    typedef std::vector<FieldEntry>::iterator iterator;                         ///< Iterator over storage vector
    typedef std::vector<FieldEntry>::const_iterator citerator;                  ///< Constant iterator over storage vector
    iterator  begin() {return Fields.begin();};                                 ///< Iterator to the begin of storage vector
    iterator  end()   {return Fields.end();};                                   ///< Iterator to the end of storage vector
    citerator cbegin() const {return Fields.cbegin();};                         ///< Constant iterator to the begin of storage vector
    citerator cend()   const {return Fields.cend();};                           ///< Constant iterator to the end of storage vector
    iterator erase(iterator it)                                                 ///< Erase a single record pointed by iterator it
    {
        it = Fields.erase(it);
        return it;
    }
    FieldEntry& front(void)                                                     ///< Reference to the first FieldEntry.
    {
        return Fields.front();
    };
    const FieldEntry& front(void) const                                         ///< Constant reference to the first FieldEntry.
    {
        return Fields.front();
    };
    int flag;                                                                   ///< Interface flag: 1 if node is near the interface, 2 if it is in the interface and 0 if it is in the bulk

 protected:
 private:
    std::vector<FieldEntry> Fields;                                             ///< Fields storage vector.
};

/***************************************************************/
inline int Node::finalize(void)
{
    double total = 0.0;
    for(auto i = Fields.begin(); i != Fields.end();)
    {
        int erase = 0;
        if (i->value >= 1.0 - DBL_EPSILON) i->value = 1.0;
        if (i->value <= 0.0 + DBL_EPSILON)
        {
            i->value = 0.0; erase = 1;
        }

        total += i->value;
        if (erase)
        {
            i = Fields.erase(i);
        }
        else
        {
            ++i;
        }
    }

    if(Fields.size() > 1)
    {
        double norm = 1.0/total;
        for(auto i = Fields.begin(); i < Fields.end(); ++i)
        {
            i->value *= norm;
            i->value2 = 0.0;
        }
        flag = 2;
    }
    else
    {
        Fields.front().value = 1.0;
        Fields.front().value2 = 0.0;
        flag = 0;
    }
    return flag;
}

/***************************************************************/
inline void Node::set(const int n, const double value)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->value = value;
        return;
    }

    if(Fields.size())
    {
        flag = 2;
    }

    FieldEntry NewEntry;
    NewEntry.index  = n;
    NewEntry.indexB = 0;
    NewEntry.value  = value;
    NewEntry.value2 = 0.0;
    Fields.push_back(NewEntry);
}

inline double Node::get(const int n) const
{
    for(auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if(i->index == n) return i->value;
    }
    return 0.0;
}

inline void Node::add(const int n, const double value)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->value += value;
        return;
    }
    if(Fields.size())
    {
        flag = 2;
    }
    FieldEntry NewEntry;
    NewEntry.index  = n;
    NewEntry.indexB = 0;
    NewEntry.value  = value;
    NewEntry.value2 = 0.0;
    Fields.push_back(NewEntry);
}

inline void Node::set2(const int n, const double value2)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->value2 = value2;
        return;
    }

    FieldEntry NewEntry;
    NewEntry.index  = n;
    NewEntry.indexB = 0;
    NewEntry.value  = 0.0;
    NewEntry.value2 = value2;
    Fields.push_back(NewEntry);
}

inline double Node::get2(const int n) const
{
    for(auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if(i->index == n) return i->value2;
    }
    return 0.0;
}

inline void Node::add2(const int n, const double value2)
{
    for(auto i = Fields.begin(); i < Fields.end(); ++i)
    if(i->index == n)
    {
        i->value2 += value2;
        return;
    }

    FieldEntry NewEntry;
    NewEntry.index  = n;
    NewEntry.indexB = 0;
    NewEntry.value  = 0.0;
    NewEntry.value2 = value2;
    Fields.push_back(NewEntry);
}
/***************************************************************/
inline double Node::operator[](const int n) const
{
    for(auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        if(i->index == n) return i->value;
    }
    return 0.0;
}

inline void Node::set(const int n, const int m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->value = value;
        return;
    }

    FieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value  = value;
    NewEntry.value2 = 0.0;

    Fields.push_back(NewEntry);

}

inline void Node::add(const int n, const int m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->value += value;
        return;
    }

    FieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value  = value;
    NewEntry.value2 = 0.0;

    Fields.push_back(NewEntry);
}

/***************************************************************/
inline double Node::get(const int n, const int m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        return i->value;
    }
    return 0.0;
}

inline void Node::set2(const int n, const int m, const double value2)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->value2 = value2;
        return;
    }

    FieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value  = 0.0;
    NewEntry.value2 = value2;

    Fields.push_back(NewEntry);

}

inline void Node::add2(const int n, const int m, const double value2)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        i->value2 += value2;
        return;
    }

    FieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value  = 0.0;
    NewEntry.value2 = value2;

    Fields.push_back(NewEntry);
}

/***************************************************************/
inline double Node::get2(const int n, const int m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if (i->indexA == n && i->indexB == m)
    {
        return i->value2;
    }
    return 0.0;
}

inline Node Node::operator+(const Node& n) const
{
    Node result = n;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add(i->indexA, i->indexB, i->value);
    }
    result.flag = std::max(flag, n.flag);
    return result;
}

inline Node Node::operator-(const Node& n) const
{
    Node result = n;

    for (auto i = result.begin(); i < result.end(); ++i)
    {
        i->value = -i->value;
    }

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add(i->indexA, i->indexB, i->value);
    }
    result.flag = std::max(flag, n.flag);
    return result;
}

inline Node Node::operator*(const double n) const
{
    Node result = *this;

    for(auto i = result.begin(); i < result.end(); ++i)
    {
        i->value *= n;
    }

    return result;
}

inline Node& Node::operator=(const Node& n)
{
    Fields = n.Fields;
    flag = n.flag;

    return *this;
}

inline void Node::set_sym(const int n, const int m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value = value;
        return;
    }

    FieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value  = value;
    NewEntry.value2 = 0.0;

    Fields.push_back(NewEntry);
}

inline void Node::set_sym2(const int n, const int m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value2 = value;
        return;
    }

    FieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value  = 0.0;
    NewEntry.value2 = value;

    Fields.push_back(NewEntry);
}

inline void Node::add_sym(const int n, const int m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        i->value += value;
        return;
    }

    FieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value  = value;
    NewEntry.value2 = 0.0;

    Fields.push_back(NewEntry);
}

/***************************************************************/
inline double Node::get_sym(const int n, const int m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return i->value;
    }
    return 0.0;
}

inline double Node::get_sym2(const int n, const int m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if ((i->indexA == n and i->indexB == m) or
        (i->indexA == m and i->indexB == n) )
    {
        return i->value2;
    }
    return 0.0;
}

inline void Node::set_asym(const int n, const int m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->value = value;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->value = -value;
        return;
    }

    FieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value  = value;
    NewEntry.value2 = 0.0;

    Fields.push_back(NewEntry);
}

inline void Node::add_asym(const int n, const int m, const double value)
{
    for (auto i = Fields.begin(); i < Fields.end(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        i->value += value;
        return;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        i->value -= value;
        return;
    }

    FieldEntry NewEntry;
    NewEntry.indexA = n;
    NewEntry.indexB = m;
    NewEntry.value  = value;
    NewEntry.value2 = 0.0;

    Fields.push_back(NewEntry);
}

/***************************************************************/
inline double Node::get_asym(const int n, const int m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        return i->value;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        return -i->value;
    }
    return 0.0;
}

/***************************************************************/
inline double Node::get_asym2(const int n, const int m) const
{
    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if (i->indexA == n and i->indexB == m)
    {
        return i->value2;
    }
    else if(i->indexA == m and i->indexB == n)
    {
        return -i->value2;
    }
    return 0.0;
}

inline Node Node::add(const Node& n) const
{
    Node result = n;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add(i->indexA, i->indexB, i->value);
    }
    result.flag = std::max(flag, n.flag);
    return result;
}

inline Node Node::add_sym(const Node& n) const
{
    Node result = n;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_sym(i->indexA, i->indexB, i->value);
    }
    return result;
}

inline Node Node::add_asym(const Node& n) const
{
    Node result = n;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        result.add_asym(i->indexA, i->indexB, i->value);
    }
    return result;
}

inline Node Node::add_exist(const Node& n) const
{
    Node result;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        double locVal = n.get(i->indexA, i->indexB);
        if(locVal != 0.0 && i->value != 0.0)
        {
            result.set(i->indexA, i->indexB, i->value + locVal);
        }
    }
    result.flag = std::max(flag, n.flag);

    return result;
}

inline Node Node::add_sym_exist(const Node& n) const
{
    Node result;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        double locVal = n.get_sym(i->indexA, i->indexB);
        if(locVal != 0.0 && i->value != 0.0)
        {
            result.set_sym(i->indexA, i->indexB, i->value + locVal);
        }
    }

    return result;
}

inline Node Node::add_asym_exist(const Node& n) const
{
    Node result;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    {
        double locVal = n.get_asym(i->indexA, i->indexB);
        if(locVal != 0.0 && i->value != 0.0)
        {
            result.set_asym(i->indexA, i->indexB, i->value + locVal);
        }
    }

    return result;
}

inline FieldEntry Node::get_max(void) const
{
    FieldEntry returnFieldEntry;

    returnFieldEntry.value = 0.0;

    for (auto i = Fields.cbegin(); i < Fields.cend(); ++i)
    if (fabs(i->value) > returnFieldEntry.value)
    {
        returnFieldEntry.index  = i->index;
        returnFieldEntry.indexB = i->indexB;
        returnFieldEntry.value  = i->value;
        returnFieldEntry.value2 = i->value2;
    }
    return returnFieldEntry;
}

}//namespace opensim
#endif