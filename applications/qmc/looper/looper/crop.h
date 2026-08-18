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

#ifndef LOOPER_CROP_H
#define LOOPER_CROP_H

#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>

namespace looper {

using boost::enable_if;
using boost::disable_if;
using boost::is_arithmetic;

//
// function crop_0, crop_01
//

template<typename T>
T crop_0(T x, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return (x > T(0)) ? x : T(0); }

template<typename T>
T crop_0(T const& x, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return (x > T(0)) ? x : T(0); }

template<typename T>
T crop_01(T x, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return (x < T(1)) ? crop_0(x) : T(1); }

template<typename T>
T crop_01(T const& x, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return (x < T(1)) ? crop_0(x) : T(1); }

} // end namespace looper

#endif // LOOPER_CROP_H
