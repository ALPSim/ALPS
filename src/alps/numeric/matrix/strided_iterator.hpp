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

#ifndef ALPS_STRIDED_ITERATOR_HPP
#define ALPS_STRIDED_ITERATOR_HPP

#include <boost/iterator/iterator_facade.hpp>
#include <boost/static_assert.hpp>
#include <cassert>

namespace alps {
namespace numeric {

template <typename Matrix, typename T>
class strided_iterator : public boost::iterator_facade<
                                strided_iterator<Matrix,T>,
                                T,
                                boost::random_access_traversal_tag,
                                T&,
                                typename Matrix::difference_type
                                >
{
    public:
        typedef T value_type;

        strided_iterator()
        : ptr(0), stride(0)
        {
        }

        strided_iterator(value_type* ptr, typename Matrix::difference_type stride)
            : ptr(ptr), stride(stride)
        {
            // The value_type of the iterator must be the value_type of the matrix or const Matrix::value_type
            BOOST_STATIC_ASSERT( (boost::is_same<typename Matrix::value_type, T>::value
                                 || boost::is_same<const typename Matrix::value_type,T>::value) );
        }

        template<typename Matrix2, typename U>
        strided_iterator(strided_iterator<Matrix2,U> const& r)
            : ptr(r.ptr), stride(r.stride)
        {
            BOOST_STATIC_ASSERT( (boost::is_same<
                        typename boost::add_const<Matrix>::type,
                        typename boost::add_const<Matrix2>::type
                        >::value ));
        }

    private:
        friend class boost::iterator_core_access;
        template <typename,typename> friend class strided_iterator;

        value_type& dereference() const
        {
            assert(ptr != 0);
            return *ptr;
        }

        template <typename Matrix2,typename U>
        bool equal(strided_iterator<Matrix2,U> const& y) const
        {
            BOOST_STATIC_ASSERT( (boost::is_same<
                        typename boost::add_const<Matrix>::type,
                        typename boost::add_const<Matrix2>::type
                        >::value ));
            return (ptr == y.ptr);
        }
        void increment()
        {
            ptr+=stride;
        }
        void decrement()
        {
            ptr-=stride;
        }
        void advance(typename Matrix::difference_type n)
        {
            ptr += n*stride;
        }

        template <typename U>
        typename Matrix::difference_type distance_to(strided_iterator<Matrix,U> const& z) const
        {
            assert( stride != 0 );
            assert( (z.ptr - ptr) % stride == 0 );
            return (z.ptr - ptr)/stride;
        }

        value_type* ptr;
        typename Matrix::difference_type stride;
};

} // end namespace numeric
} // end namespace alps

#endif //ALPS_STRIDED_ITERATOR_HPP
