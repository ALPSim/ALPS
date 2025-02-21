/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 * 
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_BLOCK_MATRIX_GROUPED_SYMMETRY_H
#define MAQUIS_DMRG_BLOCK_MATRIX_GROUPED_SYMMETRY_H

#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/block_matrix/indexing.h"

//// TRAITS

template <class SymmGroup>
struct grouped_symmetry;

template <>
struct grouped_symmetry<TrivialGroup> {
    typedef TrivialGroup type;
};

template <>
struct grouped_symmetry<U1> {
    typedef TwoU1 type;
};

//// GROUPPING FUNCTIONS

inline TrivialGroup::charge group(TrivialGroup::charge c1, TrivialGroup::charge c2)
{
    return TrivialGroup::IdentityCharge;
}

inline TwoU1::charge group(U1::charge c1, U1::charge c2)
{
    TwoU1::charge R;
    R[0] = c1; R[1] = c2;
    return R;
}

template<class SymmGroup>
Index<typename grouped_symmetry<SymmGroup>::type> group(Index<SymmGroup> const & i1,
                                             Index<SymmGroup> const & i2)
{
    typedef typename grouped_symmetry<SymmGroup>::type OutSymm;
    
    Index<OutSymm> ret;
    for (typename Index<SymmGroup>::const_iterator it1 = i1.begin(); it1 != i1.end(); ++it1)
        for (typename Index<SymmGroup>::const_iterator it2 = i2.begin(); it2 != i2.end(); ++it2)
        {
            ret.insert(std::make_pair(group(it1->first, it2->first), it1->second*it2->second));
        }
    return ret;
}


#endif
