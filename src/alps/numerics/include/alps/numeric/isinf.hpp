/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2011 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_NUMERIC_ISINF_HPP
#define ALPS_NUMERIC_ISINF_HPP

#include <alps/config.h>
#include <cmath>

namespace alps { namespace numeric {

#ifdef isinf
#undef isinf
#endif

#if defined( BOOST_MSVC)
  template <class T>
  bool isinf(T x) { return !_finite(x) && !_isnan(x);}
#elif (defined(__INTEL_COMPILER) && __cplusplus <= 199711L) || defined(_CRAYC) || defined(__FCC_VERSION)
  using ::isinf;
#else
  using std::isinf;
#endif

} } // end namespace

#endif // ALPS_NUMERIC_ISINF_HPP
