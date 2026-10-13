/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Matthias Troyer <troyer@comp-phys.org>,
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: obsvalue.h 3435 2009-11-28 14:45:38Z troyer $ */

#ifndef ALPS_TYPE_TRAITS_HAS_VALUE_TYPE_H
#define ALPS_TYPE_TRAITS_HAS_VALUE_TYPE_H

#include <boost/mpl/has_xxx.hpp>

// maybe we can automate this by checking for the existence of a value_type member

namespace alps {

BOOST_MPL_HAS_XXX_TRAIT_DEF(value_type)

} // end namespace alps

#endif // ALPS_TYPE_TRAITS_HAS_VALUE_TYPE_H
