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

#ifndef ALPS_NUMERIC_MATRIX_OPERATORS_MULTIPLY_SCALAR_HPP
#define ALPS_NUMERIC_MATRIX_OPERATORS_MULTIPLY_SCALAR_HPP

namespace alps {
namespace numeric {

template <typename T1, typename T2, typename Tag1>
typename multiply_return_type_helper<T1,T2>::type multiply(T1 const& t1, T2 const& t2, Tag1, tag::scalar)
{
    typename multiply_return_type_helper<T1,T2>::type r(t1);
    r *= t2;
    return r;
}

template <typename T1, typename T2, typename Tag2>
typename multiply_return_type_helper<T1,T2>::type multiply(T1 const& t1, T2 const& t2, tag::scalar, Tag2)
{
    typename multiply_return_type_helper<T1,T2>::type r(t2);
    r *= t1;
    return r;
}

} // end namespace numeric
} // end namespace alps

#endif // ALPS_NUMERIC_MATRIX_OPERATORS_MULTIPLY_SCALAR_HPP
