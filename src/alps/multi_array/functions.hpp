/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2012 by Ilia Zintchenko <iliazin@gmail.com>                       *
 *                       Jan Gukelberger                                           *
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

#ifndef ALPS_MULTI_ARRAY_FUNCTIONS_HPP
#define ALPS_MULTI_ARRAY_FUNCTIONS_HPP

#include <alps/multi_array/multi_array.hpp>
#include <alps/numeric/special_functions.hpp>

namespace alps {

#define ALPS_IMPLEMENT_FUNCTION(FUN)                                                                                            \
    template <class T, std::size_t D, class Allocator> multi_array<T,D,Allocator> FUN (multi_array<T, D, Allocator> arg) {      \
        using std:: FUN;                                                                                                        \
        std::transform(arg.data(), arg.data() + arg.num_elements(), arg.data(), std::ptr_fun<T, T>(FUN));                       \
        return arg;                                                                                                             \
    }

ALPS_IMPLEMENT_FUNCTION(sin)
ALPS_IMPLEMENT_FUNCTION(cos)
ALPS_IMPLEMENT_FUNCTION(tan)
ALPS_IMPLEMENT_FUNCTION(sinh)
ALPS_IMPLEMENT_FUNCTION(cosh)
ALPS_IMPLEMENT_FUNCTION(tanh)
ALPS_IMPLEMENT_FUNCTION(asin)
ALPS_IMPLEMENT_FUNCTION(acos)
ALPS_IMPLEMENT_FUNCTION(atan)
ALPS_IMPLEMENT_FUNCTION(abs)
ALPS_IMPLEMENT_FUNCTION(sqrt)
ALPS_IMPLEMENT_FUNCTION(exp)
ALPS_IMPLEMENT_FUNCTION(log)
ALPS_IMPLEMENT_FUNCTION(fabs)

#undef ALPS_IMPLEMENT_FUNCTION

#define ALPS_IMPLEMENT_FUNCTION(FUN)                                                                                            \
    template <class T, std::size_t D, class Allocator> multi_array<T, D, Allocator> FUN (multi_array<T, D, Allocator> arg) {    \
        using alps::numeric:: FUN ;                                                                                             \
        std::transform(arg.data(), arg.data() + arg.num_elements(), arg.data(), std::ptr_fun<T, T>(FUN));                       \
        return arg;                                                                                                             \
    }

ALPS_IMPLEMENT_FUNCTION(sq)
ALPS_IMPLEMENT_FUNCTION(cb)
ALPS_IMPLEMENT_FUNCTION(cbrt)

#undef ALPS_IMPLEMENT_FUNCTION

    template <class T1, class T2, std::size_t D, class Allocator>
    multi_array<T1,D,Allocator> pow(multi_array<T1,D,Allocator> a, T2 s)
  {
    std::pointer_to_binary_function <T1,T2,T1> PowObject (std::ptr_fun<T1,T2,T1>(std::pow));
    std::transform(a.data(),a.data()+a.num_elements(),a.data(),std::bind2nd(PowObject,s));
    return a;
  }

  template <class T, std::size_t D, class Allocator>
  T sum(multi_array<T,D,Allocator>& a)
  {
    return std::accumulate(a.data(),a.data()+a.num_elements(),0.,std::plus<T>());
  }

}//namespace alps

#endif // ALPS_MULTI_ARRAY_FUNCTIONS_HPP

