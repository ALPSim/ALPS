/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Andreas Hehn <hehn@phys.ethz.ch>                   *
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
#ifndef ALPS_NUMERIC_MATRIX_GEMV_HPP
#define ALPS_NUMERIC_MATRIX_GEMV_HPP

#include <alps/numeric/matrix/detail/debug_output.hpp>
#include <alps/numeric/matrix/is_blas_dispatchable.hpp>
#include <boost/numeric/bindings/blas/level2/gemv.hpp>
#include <algorithm>
#include <cassert>

namespace alps {
namespace numeric {

template <typename Matrix, typename Vector, typename Vector2>
void gemv(Matrix const& m, Vector const& x, Vector2& y, boost::mpl::false_)
{
    typedef typename Matrix::size_type size_type;
    using std::fill;
    assert(num_cols(m) > 0);
    fill(y.begin(), y.end(), typename Vector2::value_type(0));
    for(size_type j=0; j < num_cols(m); ++j)
    for(size_type i=0; i < num_rows(m); ++i)
        y[i] += m(i,j) * x[j];
}

template <typename Matrix, typename Vector, typename Vector2>
void gemv(Matrix const& m, Vector const& x, Vector2& y, boost::mpl::true_)
{
    ALPS_NUMERIC_MATRIX_DEBUG_OUTPUT( "using blas gemv for " << typeid(m).name() << " " << typeid(x).name() << " -> " << typeid(y).name() );
    typedef typename Matrix::value_type value_type;
    boost::numeric::bindings::blas::gemv(value_type(1), m, x, value_type(0), y);
}

template <typename Matrix, typename Vector, typename Vector2>
void gemv(Matrix const& m, Vector const& x, Vector2& y)
{
    assert(num_rows(m) == y.size());
    assert(num_cols(m) == x.size());
    // TODO test also Vector2
    gemv(m,x,y,is_blas_dispatchable<Matrix,Vector>());
}

} // end namespace numeric
} // end namespace alps

#endif // ALPS_NUMERIC_MATRIX_GEMV_HPP
