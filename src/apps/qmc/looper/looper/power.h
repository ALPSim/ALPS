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

#ifndef LOOPER_POWER_H
#define LOOPER_POWER_H

#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <complex>

namespace looper {

using boost::enable_if;
using boost::disable_if;
using boost::is_arithmetic;

//
// function power2, power3, power4
//

#ifndef BOOST_NO_SFINAE

template<typename T>
T power2(T t, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return t * t; }

template<typename T>
T power2(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return t * t; }

template<typename T>
T power2(std::complex<T> const& t)
{ return power2(real(t)) + power2(imag(t)); }

template<typename T>
T power3(T t, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return t * t * t; }

template<typename T>
T power3(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return t * t * t; }

template<typename T>
T power4(T t, typename enable_if<is_arithmetic<T> >::type* = 0)
{ return power2(power2(t)); }

template<typename T>
T power4(T const& t, typename disable_if<is_arithmetic<T> >::type* = 0)
{ return power2(power2(t)); }

template<typename T>
T power4(std::complex<T> const& t)
{ return power2(power2(t)); }

#else

template<typename T>
T power2(T const& t) { return t * t; }

template<typename T>
T power3(T const& t) { return t * t * t; }

template<typename T>
T power4(T const& t) { return power2(power2(t)); }

#endif

} // end namespace looper

#endif // LOOPER_POWER_H
