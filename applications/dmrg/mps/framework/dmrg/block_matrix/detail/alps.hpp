/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MAQUIS_BLOCK_MATRIX_DEATAIL_ALPS_HPP
#define MAQUIS_BLOCK_MATRIX_DEATAIL_ALPS_HPP

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>
#include <alps/numeric/diagonal_matrix.hpp>
#include "dmrg/block_matrix/detail/alps_detail.hpp"
#include "dmrg/block_matrix/detail/one_matrix.hpp"
#include "utils/traits.hpp"

namespace maquis { namespace traits {

    template <typename T, typename MemoryBlock> 
    struct transpose_view< alps::numeric::matrix<T, MemoryBlock> > { typedef alps::numeric::transpose_view<alps::numeric::matrix<T, MemoryBlock> > type; }; 

    template <typename T> 
    struct transpose_view< alps::numeric::diagonal_matrix<T> > { typedef alps::numeric::diagonal_matrix<T> type; };

} }

namespace alps { namespace numeric {

    template<class Matrix> struct associated_one_matrix { };
    template<class Matrix> struct associated_dense_matrix { };

    template<typename T, typename MemoryBlock>
    struct associated_one_matrix<alps::numeric::matrix<T, MemoryBlock> > { typedef maquis::dmrg::one_matrix<T> type; };

    template<typename T, class MemoryBlock>
    struct associated_dense_matrix<alps::numeric::matrix<T, MemoryBlock> > { typedef alps::numeric::matrix<T, MemoryBlock> type; };

    template<typename T>
    struct associated_dense_matrix<maquis::dmrg::one_matrix<T> > { typedef alps::numeric::matrix<T> type; };

} }

#endif
