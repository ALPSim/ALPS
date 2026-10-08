/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/alea/nan.h>
#include <limits>

namespace alps {

double nan() { return std::numeric_limits<double>::signaling_NaN(); }

double inf() { return std::numeric_limits<double>::infinity(); }

double ninf() { return -std::numeric_limits<double>::infinity(); }

} // namespace alps
