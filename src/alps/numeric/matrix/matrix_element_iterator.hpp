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

#ifndef ALPS_MATRIX_ELEMENT_ITERATOR_HPP
#define ALPS_MATRIX_ELEMENT_ITERATOR_HPP

#include <boost/iterator/iterator_facade.hpp>
#include <boost/static_assert.hpp>
#include <iostream>

namespace alps {
namespace numeric {
template <typename Matrix, typename T>
class matrix_element_iterator : public boost::iterator_facade<
                                matrix_element_iterator<Matrix,T>,
                                T,
                                boost::random_access_traversal_tag,
                                T&,
                                typename Matrix::difference_type
                                >
{
    public:
        typedef T value_type;

        matrix_element_iterator(Matrix* m, typename Matrix::size_type row_p, typename Matrix::size_type col_p)
            : m(m), i(row_p), j(col_p)
        {
            // The value_type of the iterator must be the value_type of the matrix or const Matrix::value_type
            BOOST_STATIC_ASSERT( (boost::is_same<typename Matrix::value_type, T>::value
                                 || boost::is_same<const typename Matrix::value_type,T>::value) );
#ifndef ALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING
            std::cerr<<"WARNING: matrix_element_iterators are very slow!"<<std::endl;
            std::cerr<<"You should use strided_iterators (eg. row_iterator) instead, unless you really don't care."<<std::endl;
            std::cerr<<"To disable this warning compile with -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING ."<<std::endl;
#endif //ALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING
        }

        template<typename Matrix2, typename U>
        matrix_element_iterator(matrix_element_iterator<Matrix2,U> const& r)
            : m(r.m), i(r.i), j(r.j)
        {
            BOOST_STATIC_ASSERT( (boost::is_same<
                        typename boost::add_const<Matrix>::type,
                        typename boost::add_const<Matrix2>::type
                        >::value ));
        }

    private:
        friend class boost::iterator_core_access;
        template <typename,typename> friend class matrix_element_iterator;

        value_type& dereference() const
        { return m->operator()(i,j); }

        template <typename Matrix2,typename U>
        bool equal(matrix_element_iterator<Matrix2,U> const& y) const
        {
            BOOST_STATIC_ASSERT( (boost::is_same<
                        typename boost::add_const<Matrix>::type,
                        typename boost::add_const<Matrix2>::type
                        >::value ));
            if( i == y.i && j == y.j )
                return true;
            else
                return false;
        }
        void increment()
        {
            ++i;
            if( i >= static_cast<typename Matrix::difference_type>(num_rows(*m)) )
            {
                i=0;
                ++j;
            }
        }
        void decrement()
        {
            --i;
            if( i < 0)
            {
                i = num_rows(*m)-1;
                --j;
            }
        }
        void advance(typename Matrix::difference_type n)
        {
            j += n / num_rows(*m);
            i += n % num_rows(*m);
        }

        template <typename U>
        typename Matrix::difference_type distance_to(matrix_element_iterator<Matrix,U> const& z) const
        {
            return (z.j - j)*num_rows(*m) + z.i - i;
        }

        // Matrix
        Matrix* m;
        // position row i, column j
        typename Matrix::difference_type i;
        typename Matrix::difference_type j;
};

} // end namespace numeric
} // end namespace alps

#endif //ALPS_MATRIX_ELEMENT_ITERATOR_HPP
