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

#ifndef ALPS_NUMERIC_MATRIX_DETAIL_AUTO_DEDUCE_PLUS_RETURN_TYPE_HPP
#define ALPS_NUMERIC_MATRIX_DETAIL_AUTO_DEDUCE_PLUS_RETURN_TYPE_HPP

#include <boost/mpl/if.hpp>

namespace alps {
namespace numeric {
namespace detail {

template <typename T1, typename T2>
struct auto_deduce_plus_return_type
{
    private:
        typedef char one;
        typedef long unsigned int two;
        static one test(T1 t) {return one();}
        static two test(T2 t) {return two();}
    public:
        typedef boost::mpl::bool_<(sizeof(test(T1()+T2())) == sizeof(one))> select_first;
        typedef typename boost::mpl::if_<select_first,T1,T2>::type type;
};

template <typename T>
struct auto_deduce_plus_return_type<T,T>
{
    typedef boost::mpl::bool_<true> select_first;
    typedef T type;
};

} // end namespace detail
} // end namespace numeric
} // end namespace alps
#endif // ALPS_NUMERIC_MATRIX_DETAIL_AUTO_DEDUCE_PLUS_RETURN_TYPE_HPP
