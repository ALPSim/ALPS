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
