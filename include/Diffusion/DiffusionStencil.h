#ifndef DIFFUSIONSTENCIL_H
#define DIFFUSIONSTENCIL_H

#include "Tools/Includes.h"

namespace opensim
{

class DiffusionStencil                                                          ///< Diffusion stencil class (uses user specified Laplacian stencil as the basis). Allows replacing the loop over Laplacian elements by the iterator which is beneficial for compact stencils 
{
 public:
    struct StencilEntry
    {
        int di;                                                                 ///< x coordinate of stencil element
        int dj;                                                                 ///< y coordinate of stencil element
        int dk;                                                                 ///< z coordinate of stencil element
        double weight;                                                          ///< weight associated with the stencil element
    };
    void Set(const double LaplacianStencil[3][3][3], double dx)                 ///< Sets the diffusion stencil using user specified Laplacian stencil
    {
        if(not StencilElements.empty())
        {
            StencilElements.clear();
        }
        
        for(int x = -1; x <= 1; ++x)
        for(int y = -1; y <= 1; ++y)
        for(int z = -1; z <= 1; ++z)
        if(x != 0 or y != 0 or z!= 0)
        if (LaplacianStencil[x+1][y+1][z+1] != 0.0)
        {
            StencilEntry locElement;
            locElement.di = x;
            locElement.dj = y;
            locElement.dk = z;
            
            locElement.weight = LaplacianStencil[x+1][y+1][z+1]/(dx*dx);
            StencilElements.push_back(locElement);        
        }    
    };
    
    unsigned int size() const {return StencilElements.size();};                 ///< Returns the size of stencil.
    typedef std::vector<StencilEntry>::iterator iterator;                       ///< Iterator over stencil
    typedef std::vector<StencilEntry>::const_iterator citerator;                ///< Constant iterator over stencil
    iterator  begin() {return StencilElements.begin();};                        ///< Iterator to the begin of sencil
    iterator  end()   {return StencilElements.end();};                          ///< Iterator to the end of stencil
    citerator cbegin() const {return StencilElements.cbegin();};                ///< Constant iterator to the begin of stencil
    citerator cend()   const {return StencilElements.cend();};                  ///< Constant iterator to the end of stencil
    iterator erase(iterator it) {it = StencilElements.erase(it); return it;}    ///< Erases entry
 protected:
 private:
 std::vector< StencilEntry > StencilElements;                                   ///< Stencil storage
};
}
#endif//DIFFUSIONSTENCIL_H