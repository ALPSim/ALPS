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

#ifndef OSIRIS_BOOST_ARRAY_HPP
#define OSIRIS_BOOST_ARRAY_HPP

#include <alps/config.h>
#include <alps/osiris/dump.h>

#include <boost/array.hpp>

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
#endif

template<class T, std::size_t N> 
inline alps::IDump& operator>>(alps::IDump& dump, boost::array<T, N>& x)
{
  dump.read_array(N,&(x[0]));
  return dump;
}

template<class T, std::size_t N> 
inline alps::ODump& operator<<(alps::ODump& dump, const boost::array<T,N>& x)
{
  dump.write_array(N,&(x[0]));
  return dump;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace alps
#endif

#endif // OSIRIS_BOOST_ARRAY_HPP
