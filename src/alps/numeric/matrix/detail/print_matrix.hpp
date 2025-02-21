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
#ifndef ALPS_NUMERIC_MATRIX_DETAIL_PRINT_MATRIX_HPP
#define ALPS_NUMERIC_MATRIX_DETAIL_PRINT_MATRIX_HPP

#include <ostream>


namespace alps {
namespace numeric {
namespace detail {

template <typename Matrix>
void print_matrix(std::ostream& os, Matrix const& m)
{
    os << "[";
    for(typename Matrix::size_type i=0; i < num_rows(m); ++i)
    {
        os << "[ ";
        if(num_cols(m) > 0)
        {
            for(typename Matrix::size_type j=0; j < num_cols(m)-1; ++j)
                os << m(i,j) << ", ";
            os << m(i,num_cols(m)-1);
        }
        os << "]";
        if(i+1 < num_rows(m))
            os << "," << std::endl;
    }
    os << "]" << std::endl;
}

} // end namespace detail
} // end namespace numeric
} // end namespace alps

#endif // ALPS_NUMERIC_MATRIX_DETAIL_PRINT_MATRIX_HPP
