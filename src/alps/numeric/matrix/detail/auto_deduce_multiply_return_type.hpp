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

#ifndef ALPS_NUMERIC_MATRIX_DETAIL_AUTO_DEDUCE_MULTIPLY_RETURN_TYPE_HPP
#define ALPS_NUMERIC_MATRIX_DETAIL_AUTO_DEDUCE_MULTIPLY_RETURN_TYPE_HPP

#include <boost/mpl/if.hpp>

namespace alps {
namespace numeric {
namespace detail {

template <typename T1, typename T2>
struct auto_deduce_multiply_return_type
{
    private:
        typedef char one;
        typedef long unsigned int two;
        static one test(T1 t) {return one();}
        static two test(T2 t) {return two();}
    public:
        typedef boost::mpl::bool_<(sizeof(test(T1()*T2())) == sizeof(one))> select_first;
        typedef typename boost::mpl::if_<select_first,T1,T2>::type type;
};

template <typename T>
struct auto_deduce_multiply_return_type<T,T>
{
    typedef boost::mpl::bool_<true> select_first;
    typedef T type;
};

} // end namespace detail
} // end namespace numeric
} // end namespace alps
#endif // ALPS_NUMERIC_MATRIX_DETAIL_AUTO_DEDUCE_MULTIPLY_RETURN_TYPE_HPP
