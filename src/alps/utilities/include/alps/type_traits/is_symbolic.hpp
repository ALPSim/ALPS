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

/* $Id: is_symbolic.hpp 3435 2009-11-28 14:45:38Z troyer $ */

#ifndef ALPS_TYPE_TRAITS_IS_SYMBOLIC_H
#define ALPS_TYPE_TRAITS_IS_SYMBOLIC_H

#include <boost/mpl/bool.hpp>

// maybe we can automate this by checking for the existence of a value_type member

namespace alps {

template <class T>
struct is_symbolic : public boost::mpl::false_ {};


} // end namespace alps

#endif // ALPS_TYPE_TRAITS_IS_SYMBOLIC_H
