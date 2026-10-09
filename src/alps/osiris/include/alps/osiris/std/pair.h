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

#ifndef OSIRIS_STD_PAIR_H
#define OSIRIS_STD_PAIR_H

#include <alps/config.h>
#include <alps/osiris/dump.h>
#include <utility>

//=======================================================================
// pair templates
//-----------------------------------------------------------------------

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
#endif

template <class T1, class T2>
inline alps::IDump& operator>>(alps::IDump& dump, std::pair<T1,T2>& x)
{
  return dump >> x.first >> x.second;
}

template <class T1, class T2>
inline alps::ODump& operator<<(alps::ODump& dump, const std::pair<T1,T2>& x)
{
  return dump << x.first << x.second;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace alps
#endif

#endif // OSIRIS_STD_PAIR_H
