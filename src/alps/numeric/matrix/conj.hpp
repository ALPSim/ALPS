/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Andreas Hehn <hehn@phys.ethz.ch>                   *
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

#ifndef ALPS_MATRIX_CONJ_HPP
#define ALPS_MATRIX_CONJ_HPP

#include <alps/numeric/conj.hpp>

namespace alps {
namespace numeric {

/**
  * Does an conj_inplace on all elements of the matrix
  */
template <typename Matrix>
void conj_inplace(Matrix& m, tag::matrix)
{
    // TODO discuss conj() for matrix may be misleading:
    //      elementwise conj() <-> adjoint()
    //
    typedef typename Matrix::col_element_iterator col_element_iterator;
    for(typename Matrix::size_type j=0; j < num_cols(m); ++j)
       for(std::pair<col_element_iterator,col_element_iterator> range = col(m,j); range.first != range.second; ++range.first)
           conj_inplace(*range.first);
}

template <typename Vector>
void conj_inplace(Vector& t, tag::vector)
{
    typename Vector::iterator const end = t.end();
    for(typename Vector::iterator it=t.begin(); it != end; ++it)
        conj_inplace(*it);
}

/**
  * Returns a matrix or a vector of type T containing a the complex conjugates of the original object.
  * It does an element-wise conjugation.
  */
template <typename T>
typename boost::enable_if<
      boost::mpl::or_<
            boost::is_same<typename get_entity<T>::type, tag::vector>
          , boost::is_same<typename get_entity<T>::type, tag::matrix>
       >
    , T
>::type conj(T const& t)
{
    T r(t);
    conj_inplace(r);
    return r;
}

} // end namespace numeric
} // end namespace alps 

#endif //ALPS_MATRIX_CONJ_HPP
