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

using std::isinf;

} } // end namespace

#endif // ALPS_NUMERIC_ISINF_HPP
