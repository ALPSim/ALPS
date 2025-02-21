/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
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
