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

#ifndef ALPS_NUMERIC_IS_POSITIVE_HPP
#define ALPS_NUMERIC_IS_POSITIVE_HPP

#include <alps/numeric/is_zero.hpp>

namespace alps { namespace numeric {

//
// is_positive
//

template<unsigned int N, class T>
inline bool is_positive(T x,
  typename boost::enable_if<boost::is_float<T> >::type* = 0)
{ return is_nonzero<N>(x) && x > T(0); }
template<unsigned int N, class T>
inline bool is_positive(T x,
  typename boost::enable_if<boost::is_integral<T> >::type* = 0)
{ return x > T(0); }
template<unsigned int N, class T>
inline bool is_positive(const T& x,
  typename boost::disable_if<boost::is_arithmetic<T> >::type* = 0)
{ return is_nonzero<N>(x) && x > T(0); }

template<class T>
inline bool is_positive(T x,
  typename boost::enable_if<boost::is_float<T> >::type* = 0)
{ return is_nonzero(x) && x > T(0); }
template<class T>
inline bool is_positive(T x,
  typename boost::enable_if<boost::is_integral<T> >::type* = 0)
{ return x > T(0); }
template<class T>
inline bool is_positive(const T& x,
  typename boost::disable_if<boost::is_arithmetic<T> >::type* = 0)
{ return is_nonzero(x) && x > T(0); }


} } // end namespace alps::numeric

#endif // ALPS_NUMERIC_IS_POSITIVE_HPP
