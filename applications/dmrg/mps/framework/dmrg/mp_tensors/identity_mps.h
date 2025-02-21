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

#ifndef MAQUIS_DMRG_IDENTITY_MPS_H
#define MAQUIS_DMRG_IDENTITY_MPS_H

#include "dmrg/block_matrix/grouped_symmetry.h"

template <class Matrix, class InSymm>
MPS<Matrix, typename grouped_symmetry<InSymm>::type> identity_dm_mps(std::size_t L, Index<InSymm> const& phys_psi,
                                                                     std::vector<Index<typename grouped_symmetry<InSymm>::type> > const& allowed)
{
    MPOTensor<Matrix, InSymm> t(1,1);
    t.set(0, 0, identity_matrix<Matrix>(phys_psi));

    MPO<Matrix, InSymm> mpo(L, t);
    return mpo_to_smps_group(mpo, phys_psi, allowed);
}

#endif
