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
#ifndef ALPS_NUMERIC_MATRIX_OPERATORS_PLUS_MINUS_HPP
#define ALPS_NUMERIC_MATRIX_OPERATORS_PLUS_MINUS_HPP

#include <alps/numeric/matrix/entity.hpp>
#include <boost/type_traits/remove_const.hpp>

namespace alps {
namespace numeric {

template <typename T1, typename T2, typename EntityTag1, typename EntityTag2>
struct plus_minus_return_type
{
};


template <typename T1, typename T2>
struct plus_minus_return_type_helper
: plus_minus_return_type<
      typename boost::remove_const<T1>::type
    , typename boost::remove_const<T2>::type
    , typename get_entity<T1>::type
    , typename get_entity<T2>::type
> {
};

template <typename T1, typename T2, typename Tag1, typename Tag2>
typename plus_minus_return_type_helper<T1,T2>::type do_plus(T1 const& t1, T2 const& t2, Tag1, Tag2)
{
    typename plus_minus_return_type_helper<T1,T2>::type r(t1);
    r += t2;
    return r;
}

template <typename T1, typename T2, typename Tag1, typename Tag2>
typename plus_minus_return_type_helper<T1,T2>::type do_minus(T1 const& t1, T2 const& t2, Tag1, Tag2)
{
    typename plus_minus_return_type_helper<T1,T2>::type r(t1);
    r -= t2;
    return r;
}


template <typename T1, typename T2>
typename plus_minus_return_type_helper<T1,T2>::type operator + (T1 const& t1, T2 const& t2)
{
    return do_plus(t1, t2, typename get_entity<T1>::type(), typename get_entity<T2>::type());
}

template <typename T1, typename T2>
typename plus_minus_return_type_helper<T1,T2>::type operator - (T1 const& t1, T2 const& t2)
{
    return do_minus(t1, t2, typename get_entity<T1>::type(), typename get_entity<T2>::type());
}

} // end namespace numeric
} // end namespace alps

#endif // ALPS_NUMERIC_MATRIX_OPERATORS_PLUS_MINUS_HPP
