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

#ifndef OSIRIS_STD_VALARRAY_H
#define OSIRIS_STD_VALARRAY_H

#include <alps/config.h>

#include <alps/osiris/dump.h>
#include <alps/osiris/std/impl.h>
#include <valarray>

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
#endif

/// deserialize a std::valarray container
template <class T>
inline alps::IDump& operator >> (alps::IDump& dump, std::valarray<T>& x)
{
  x.resize(uint32_t(dump));
  dump.read_array(x.size(),&(x[0]));
  return dump;
}

/// serialize a std::valarray container
template <class T>
inline alps::ODump& operator << (alps::ODump& dump,
                                   const std:: valarray<T>& x)
{
  dump << uint32_t(x.size());
  dump.write_array(x.size(),&(const_cast<std::valarray<T>&>(x)[0]));
  return dump;
}          

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace alps
#endif

#endif // OSIRIS_STD_VALARRAY_H
