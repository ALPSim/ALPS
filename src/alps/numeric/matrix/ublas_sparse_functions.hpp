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
#ifndef ALPS_NUMERIC_MATRIX_UBLAS_SPARSE_FUNCTIONS_HPP
#define ALPS_NUMERIC_MATRIX_UBLAS_SPARSE_FUNCTIONS_HPP

#include <alps/numeric/matrix/entity.hpp>
#include <alps/numeric/matrix/operators/multiply.hpp>
#include <boost/numeric/ublas/fwd.hpp>

namespace alps {
namespace numeric {

template<class T, class L, class A>
struct entity< ::boost::numeric::ublas::mapped_vector_of_mapped_vector<T,L,A> >
{
    typedef tag::matrix type;
};

template <typename T, typename L, typename A, typename Vector>
typename multiply_return_type_helper<matrix<T>,Vector>::type multiply(::boost::numeric::ublas::mapped_vector_of_mapped_vector<T,L,A> const& m, Vector const& t2, tag::matrix, tag::vector)
{
    // NOTE: This function depends on some implementation details of ublas. I didn't see any other efficient way.
    using ::boost::numeric::ublas::map_std;
    BOOST_STATIC_ASSERT(( boost::is_same<A, map_std<std::size_t, map_std<std::size_t, T> > >::value ));
    BOOST_STATIC_ASSERT(( boost::is_same<L, ::boost::numeric::ublas::row_major>::value ));
    typedef A array_type;
    typedef typename map_std<std::size_t, map_std<std::size_t, T> >::const_iterator const_iterator1;
    typedef typename map_std<std::size_t, T>::const_iterator                        const_iterator2;

    // We just use the same return type as for a regular dense matrix vector multiplication
    typename multiply_return_type_helper<matrix<T>,Vector>::type r(m.size1());
    const_iterator1 const end = m.data().end();
    for(const_iterator1 col_it = m.data().begin(); col_it != end; ++col_it)
    {
        const_iterator2 const row_end = col_it->second.end();
        for(const_iterator2 row_it = col_it->second.begin(); row_it != row_end; ++row_it)
            r(col_it->first) += row_it->second * t2(row_it->first);
    }
    return r;
}

}
}

#endif // ALPS_NUMERIC_MATRIX_UBLAS_SPARSE_FUNCTIONS_HPP
