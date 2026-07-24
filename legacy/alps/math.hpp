/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_MATH_HPP
#define ALPS_MATH_HPP

#if defined(_MSC_VER) || defined(__BORLANDC__) || defined(__DMC__)
#  pragma message ("This header is deprecated. Please use the new headers in alps/numeric")
#elif defined(__GNUC__) || defined(__HP_aCC) || defined(__SUNPRO_CC) || defined(__IBMCPP__)
#  warning "This header is deprecated. Please use the new headers in alps/numeric"
#endif

#include <alps/numeric/real.hpp>
#include <alps/numeric/abs2.hpp>
#include <alps/numeric/binomial.hpp>
#include <alps/numeric/is_equal.hpp>
#include <alps/numeric/is_zero.hpp>
#include <alps/numeric/is_nonzero.hpp>
#include <alps/numeric/round.hpp>
#include <alps/numeric/is_negative.hpp>
#include <alps/numeric/is_positive.hpp>
#include <alps/numeric/double2int.hpp>

namespace alps {

using numeric::real;
using numeric::binomial;
using numeric::abs2;
using numeric::is_equal;
using numeric::is_zero;
using numeric::is_nonzero;
using numeric::round;
using numeric::is_negative;
using numeric::is_positive;
using numeric::double2int;

} // end namespace

#endif // ALPS_MATH_HPP
