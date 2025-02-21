/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013 by Andreas Hehn <hehn@phys.ethz.ch>                          *
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
#ifndef ALPS_NUMERIC_MATRIX_OPERATORS_MULTIPLY_MATRIX_HPP
#define ALPS_NUMERIC_MATRIX_OPERATORS_MULTIPLY_MATRIX_HPP

#include <alps/numeric/matrix/gemm.hpp>
#include <alps/numeric/matrix/gemv.hpp>

namespace alps {
namespace numeric {

template <typename Matrix1, typename Matrix2>
typename multiply_return_type_helper<Matrix1,Matrix2>::type multiply(Matrix1 const& m1, Matrix2 const& m2, tag::matrix, tag::matrix)
{
    typename multiply_return_type_helper<Matrix1,Matrix2>::type r(num_rows(m1),num_cols(m2));
    gemm(m1,m2,r);
    return r;
}

template <typename Matrix, typename Vector>
typename multiply_return_type_helper<Matrix,Vector>::type multiply(Matrix const& m, Vector const& v, tag::matrix, tag::vector)
{
    typename multiply_return_type_helper<Matrix,Vector>::type r(num_rows(m));
    gemv(m,v,r);
    return r;
}

} // end namespace numeric
} // end namespace alps

#endif // ALPS_NUMERIC_MATRIX_OPERATORS_MULTIPLY_MATRIX_HPP
