/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef LOOPER_DIVIDE_IF_POSITIVE_H
#define LOOPER_DIVIDE_IF_POSITIVE_H

#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>

namespace looper {

using boost::enable_if;
using boost::disable_if;
using boost::is_arithmetic;

//
// dip (divide_if_positive)
//

template<typename T, typename U>
T dip(T x, U y,
      typename enable_if<is_arithmetic<T> >::type* = 0,
      typename enable_if<is_arithmetic<U> >::type* = 0)
{ return (y > U(0)) ? (x / y) : T(0); }

template<typename T, typename U>
T dip(T x, U const& y,
      typename enable_if<is_arithmetic<T> >::type* = 0,
      typename disable_if<is_arithmetic<U> >::type* = 0)
{ return (y > U(0)) ? (x / y) : T(0); }

template<typename T, typename U>
T dip(T const& x, U y,
      typename disable_if<is_arithmetic<T> >::type* = 0,
      typename enable_if<is_arithmetic<U> >::type* = 0)
{ return (y > U(0)) ? (x / y) : T(0); }

template<typename T, typename U>
T dip(T const& x, U const& y,
      typename disable_if<is_arithmetic<T> >::type* = 0,
      typename disable_if<is_arithmetic<U> >::type* = 0)
{ return (y > U(0)) ? (x / y) : T(0); }

template<typename T, typename U>
T divide_if_positive(T x, U y,
      typename enable_if<is_arithmetic<T> >::type* = 0,
      typename enable_if<is_arithmetic<U> >::type* = 0)
{ return dip(x, y); }

template<typename T, typename U>
T divide_if_positive(T x, U const& y,
      typename enable_if<is_arithmetic<T> >::type* = 0,
      typename disable_if<is_arithmetic<U> >::type* = 0)
{ return dip(x, y); }

template<typename T, typename U>
T divide_if_positive(T const& x, U y,
      typename disable_if<is_arithmetic<T> >::type* = 0,
      typename enable_if<is_arithmetic<U> >::type* = 0)
{ return dip(x, y); }

template<typename T, typename U>
T divide_if_positive(T const& x, U const& y,
      typename disable_if<is_arithmetic<T> >::type* = 0,
      typename disable_if<is_arithmetic<U> >::type* = 0)
{ return dip(x, y); }

} // end namespace looper

#endif // LOOPER_DIVIDE_IF_POSITIVE_H
