/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2012 by Andreas Hehn <hehn@phys.ethz.ch>                          *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
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
