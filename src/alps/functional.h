/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2005 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

/// \file functional.h
/// \brief extensions to the standard functional header
///
/// This header contains mathematical function objects not present in
/// the standard or boost libraries.

#ifndef ALPS_FUNCTIONAL_H
#define ALPS_FUNCTIONAL_H

#include <alps/numeric/matrix/detail/auto_deduce_multiply_return_type.hpp>
#include <alps/numeric/conj.hpp>

namespace alps {


template <class T1, class T2>
struct conj_mult_return_type : ::alps::numeric::detail::auto_deduce_multiply_return_type<T1,T2> { };

/// \brief a function object for conj(x)*y
///
/// the version for real data types is just the same as std::multiplies
/// \param T the type of the arguments and result
template <class T1, class T2>
struct conj_mult
{
/// \brief returns x*y
  typename conj_mult_return_type<T1,T2>::type operator()(const T1& a, const T2& b) const { return a*b; }
};

/// \brief a function object for conj(x)*y
///
/// the version for complex data types is specialized
/// \param T the type of the arguments and result
template <class T1, class T2>
struct conj_mult<std::complex<T1>,T2>
{
/// \brief returns std::conj(x)*y
  typename conj_mult_return_type<std::complex<T1>,T2>::type operator()(const std::complex<T1>& a, const T2& b) const {
  return std::conj(a)*b;
}
};

}

#endif // ALPS_FUNCTIONAL_H
