/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef OSIRIS_STD_LIST_H
#define OSIRIS_STD_LIST_H

#include <alps/config.h>
#include <alps/osiris/dump.h>
#include <alps/osiris/std/impl.h>
#include <alps/osiris/std/pair.h>

#include <list>

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
#endif

/// deserialize a std::list container
template <class T, class Allocator>
inline alps::IDump& operator>>(alps::IDump& dump, std::list<T,Allocator>& x)
{
  alps::detail::loadContainer(dump,x);
  return dump;
}

/// serialize a std::list container
template <class T, class Allocator>
inline alps::ODump& operator<<(alps::ODump& dump, const std::list<T,Allocator>& x)
{
  alps::detail::saveContainer(dump,x);
  return dump;
}          

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace alps
#endif

#endif // OSIRIS_STD_LIST_H
