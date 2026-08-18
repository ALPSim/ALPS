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

#ifndef ALPS_TYPE_TRAITS_NORM_TYPE_H
#define ALPS_TYPE_TRAITS_NORM_TYPE_H

#include <alps/type_traits/real_type.hpp>

// maybe we can automate this by checking for the existence of a value_type member

namespace alps {

template <class T>
struct norm_type  : public real_type<T> {};

} // end namespace alps

#endif // ALPS_TYPE_TRAITS_NORM_TYPE_H
