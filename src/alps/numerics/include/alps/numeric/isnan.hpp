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

#ifndef ALPS_NUMERIC_ISNAN_HPP
#define ALPS_NUMERIC_ISNAN_HPP

#include <cmath>
#include <alps/config.h>

namespace alps { namespace numeric {

#ifdef isnan
#undef isnan
#endif

#if defined(BOOST_MSVC)
  template <class T>
  bool isnan(T x) { return _isnan(x);}
#elif (defined(__INTEL_COMPILER) && __cplusplus <= 199711L) || defined(_CRAYC) || defined(__FCC_VERSION)
  using ::isnan;
#else
  using std::isnan;
#endif

} } // end namespace

#endif // ALPS_NUMERIC_ISNAN_HPP
