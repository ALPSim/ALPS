/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2010 by Bela Bauer <bauerb@cnsi.ucsb.edu>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: real.hpp 3958 2010-03-05 09:24:06Z troyer $ */

#ifndef ALPS_NUMERIC_IMAG_HPP
#define ALPS_NUMERIC_IMAG_HPP

#include <algorithm>
#include <complex>
#include <vector>

namespace alps { namespace numeric {


template <class T>
inline T imag(T x) { return T(0);}

template <class T>
inline T imag(std::complex<T> x) { return std::imag(x); }

} }  // end namespace alps::numeric

#endif // ALPS_MATH_HPP
