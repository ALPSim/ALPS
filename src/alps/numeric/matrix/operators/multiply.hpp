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
#ifndef ALPS_NUMERIC_MATRIX_OPERATORS_MULTIPLY_HPP
#define ALPS_NUMERIC_MATRIX_OPERATORS_MULTIPLY_HPP

#include <alps/numeric/matrix/detail/auto_deduce_multiply_return_type.hpp>
#include <alps/numeric/matrix/entity.hpp>
#include <alps/numeric/matrix/exchange_value_type.hpp>
#include <boost/type_traits/remove_const.hpp>

namespace alps {
namespace numeric {

template <typename T1, typename T2, typename EntityTag1, typename EntityTag2>
struct multiply_return_type
{
};

template <typename T1, typename T2>
struct multiply_return_type<T1,T2,tag::matrix,tag::matrix>
{
};

template <typename T1, typename T2>
struct multiply_return_type<T1,T2,tag::matrix,tag::vector>
{
    typedef typename exchange_value_type<T2,typename detail::auto_deduce_multiply_return_type<typename T1::value_type,typename T2::value_type>::type>::type type;
};

template <typename T1, typename T2>
struct multiply_return_type<T1,T2,tag::matrix,tag::scalar>
{
    typedef T1 type;
};

template <typename T1, typename T2>
struct multiply_return_type<T1,T2,tag::scalar,tag::matrix>
: multiply_return_type<T2,T1,tag::matrix,tag::scalar>
{
};

template <typename T1, typename T2>
struct multiply_return_type<T1,T2,tag::vector,tag::scalar>
{
    typedef T1 type;
};

template <typename T1, typename T2>
struct multiply_return_type<T1,T2,tag::scalar,tag::vector>
: multiply_return_type<T2,T1,tag::vector,tag::scalar>
{
};

template <typename T1, typename T2>
struct multiply_return_type_helper
: multiply_return_type<
      typename boost::remove_const<T1>::type
    , typename boost::remove_const<T2>::type
    , typename get_entity<T1>::type
    , typename get_entity<T2>::type
> {
};


template <typename T1, typename T2>
typename multiply_return_type_helper<T1,T2>::type operator * (T1 const& t1, T2 const& t2)
{
    return multiply(t1, t2, typename get_entity<T1>::type(), typename get_entity<T2>::type());
}

} // end namespace numeric
} // end namespace alps

#endif // ALPS_NUMERIC_MATRIX_OPERATORS_MULTIPLY_HPP
