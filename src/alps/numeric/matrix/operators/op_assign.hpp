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
#ifndef ALPS_NUMERIC_OPERATORS_OP_ASSIGN_HPP
#define ALPS_NUMERIC_OPERATORS_OP_ASSIGN_HPP

#include <boost/static_assert.hpp>

namespace alps {
namespace numeric {

template <typename T>
struct not_implemented
{
    static bool const value = false;
};

template <typename T1, typename T2, typename Category1, typename Category2>
void plus_assign(T1& t1, T2 const& t2, Category1, Category2)
{
    BOOST_STATIC_ASSERT(not_implemented<T1>::value);
}

template <typename T1, typename T2, typename Category1, typename Category2>
void minus_assign(T1& t1, T2 const& t2, Category1, Category2)
{
    BOOST_STATIC_ASSERT(not_implemented<T1>::value);
}

template <typename T1, typename T2, typename Category1, typename Category2>
void multiplies_assign(T1& t1, T2 const& t2, Category1, Category2)
{
    BOOST_STATIC_ASSERT(not_implemented<T1>::value);
}

} // end namespace numeric
} // end namespace alps

#endif // ALPS_NUMERIC_OPERATORS_OP_ASSIGN_HPP
