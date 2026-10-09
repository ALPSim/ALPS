/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
