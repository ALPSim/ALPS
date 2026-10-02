/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2005 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
