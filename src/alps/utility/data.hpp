/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1999-2010 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */


#ifndef ALPS_UTILITY_DATA_HPP
#define ALPS_UTILITY_DATA_HPP

#include <alps/type_traits/element_type.hpp>

namespace alps {

/// returns a pointer to the start of storage of a container
template <class C>
inline typename element_type<C>::type* data(C& c) { return &c[0];}

/// returns a pointer to the start of storage of a container
template <class C>
inline const typename element_type<C>::type* data(const C& c) { return &const_cast<C&>(c)[0];}

} // namespace alps

#endif // ALPS_UTILITY_DATA_HPP
