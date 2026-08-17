/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2012 by Andreas Hehn <hehn@phys.ethz.ch>                          *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NUMERIC_MATRIX_TRANSPOSE_HPP
#define ALPS_NUMERIC_MATRIX_TRANSPOSE_HPP
#include <alps/numeric/matrix/transpose_view.hpp>

namespace alps {
namespace numeric {

template <typename Matrix>
inline transpose_view<Matrix> transpose(Matrix const& m) {
    return transpose_view<Matrix>(m);
}

template <typename Matrix>
void transpose_inplace(Matrix& m) {
    typedef typename Matrix::size_type size_type;
    using std::swap;
    if(num_rows(m) == num_cols(m) ) {
        for(size_type i = 0; i < num_rows(m); ++i)
            for(size_type j = i+1; j < num_cols(m); ++j)
                swap(m(i,j),m(j,i));
    } else {
        // TODO replace this code by an actual inplace implementation
        Matrix m2 = transpose(m);
        swap(m,m2);
    }
}

} // end namespace numeric
} // end namespace alps
#endif //ALPS_NUMERIC_MATRIX_TRANSPOSE_HPP
