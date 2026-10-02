/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013 by Andreas Hehn <hehn@phys.ethz.ch>                          *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
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
