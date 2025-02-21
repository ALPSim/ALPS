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

#ifndef ALPS_MULTI_ARRAY_OPERATORS_HPP
#define ALPS_MULTI_ARRAY_OPERATORS_HPP

#include <alps/multi_array/multi_array.hpp>

namespace alps{

  template <class T, std::size_t D, class Allocator>
  multi_array<T,D,Allocator> operator-(multi_array<T,D,Allocator> arg)
  {
    std::transform(arg.data(), arg.data() + arg.num_elements(), arg.data(), std::negate<T>());
    return arg;
  }

  template <class T, std::size_t D, class Allocator>
  multi_array<T,D,Allocator> operator+(multi_array<T,D,Allocator> a, const multi_array<T,D,Allocator>& b)
  {
    a += b;
    return a;
  }

  template <class T, std::size_t D, class Allocator>
  multi_array<T,D,Allocator> operator-(multi_array<T,D,Allocator> a, const multi_array<T,D,Allocator>& b)
  {
    a -= b;
    return a;
  }

  template <class T, std::size_t D, class Allocator>
  multi_array<T,D,Allocator> operator*(multi_array<T,D,Allocator> a, const multi_array<T,D,Allocator>& b)
  {
    a *= b;
    return a;
  }

  template <class T, std::size_t D, class Allocator>
  multi_array<T,D,Allocator> operator/(multi_array<T,D,Allocator> a, const multi_array<T,D,Allocator>& b)
  {
    a /= b;
    return a;
  }

  template <class T1, class T2, std::size_t D, class Allocator>
  multi_array<T1,D,Allocator> operator+(const T2& b, multi_array<T1,D,Allocator> a)
  {
    a += T1(b);
    return a;
  }
  template <class T1, class T2, std::size_t D, class Allocator>
  multi_array<T1,D,Allocator> operator+(multi_array<T1,D,Allocator> a, const T2& b)
  {
    a += T1(b);
    return a;
  }

  template <class T1, class T2, std::size_t D, class Allocator>
  multi_array<T1,D,Allocator> operator-(const T2& b, multi_array<T1,D,Allocator> a)
  {
    a -= T1(b);
    return a;
  }
  template <class T1, class T2, std::size_t D, class Allocator>
  multi_array<T1,D,Allocator> operator-(multi_array<T1,D,Allocator> a, const T2& b)
  {
    a -= T1(b);
    return -a;
  }

  template <class T1, class T2, std::size_t D, class Allocator>
  multi_array<T1,D,Allocator> operator*(const T2& b, multi_array<T1,D,Allocator> a)
  {
    a *= T1(b);
    return a;
  }

  template <class T1, class T2, std::size_t D, class Allocator>
  multi_array<T1,D,Allocator> operator*(multi_array<T1,D,Allocator> a, const T2& b)
  {
    a *= T1(b);
    return a;
  }

  template <class T1, class T2, std::size_t D, class Allocator>
  multi_array<T1,D,Allocator> operator/(multi_array<T1,D,Allocator> a, const T2& b)
  {
    a /= T1(b);
    return a;
  }
  template <class T1, class T2, std::size_t D, class Allocator>
  multi_array<T2,D,Allocator> operator/(T2 const & scalar, multi_array<T1,D,Allocator> array)
  {
    std::transform(array.data(), array.data() + array.num_elements(), array.data(), std::bind1st(std::divides<T1>(),T1(scalar)));
    return array;
  }

}//namespace alps

#endif // ALPS_MULTI_ARRAY_OPERATORS_HPP
