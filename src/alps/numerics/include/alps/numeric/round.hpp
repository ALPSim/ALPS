/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_NUMERIC_ROUND_HPP
#define ALPS_NUMERIC_ROUND_HPP

#include <alps/numeric/is_zero.hpp>
#include <complex>

namespace alps { namespace numeric {

//
// round
//

/// \brief rounds a floating point value to be exactly zero if it is nearly zero
///
/// the function is specialized for floating point and complex types
/// and does nothing for other types
/// \return 0. if the floating point value of the argument is less
/// than epsilon (1e-50 by default), and the argument itself otherwise
template<unsigned int N, class T>
inline T round(T x,
  typename boost::enable_if<boost::is_arithmetic<T> >::type* = 0)
{ return is_zero<N>(x) ? T(0) : x; }

namespace detail {

template<unsigned int N, class T>
struct round_helper
{
  static T result(const T& x) { return x; }
};
template<unsigned int N, class T>
struct round_helper<N, std::complex<T> >
{
  static std::complex<T> result(const std::complex<T>& x)
  { return std::complex<double>(round<N>(x.real()), round<N>(x.imag())); }
};

} // end namespace detail

template<unsigned int N, class T>
inline T round(const T& x,
  typename boost::disable_if<boost::is_arithmetic<T> >::type* = 0)
{ return detail::round_helper<N, T>::result(x); }

template<class T>
inline T round(T x,
  typename boost::enable_if<boost::is_arithmetic<T> >::type* = 0)
{ return is_zero(x) ? T(0) : x; }
template<class T>
inline T round(const T& x,
  typename boost::disable_if<boost::is_arithmetic<T> >::type* = 0)
{ return x; }
template<class T>
inline std::complex<T> round(const std::complex<T>& x)
{ return std::complex<T>(round(x.real()), round(x.imag())); }

} } // end namespace

#endif // ALPS_NUMERIC_ROUND_HPP
