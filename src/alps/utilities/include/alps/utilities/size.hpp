/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@comp-phys.org>,
*                            Beat Ammon <ammon@ginnan.issp.u-tokyo.ac.jp>,
*                            Andreas Laeuchli <laeuchli@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_UTILITY_SIZE_HPP
#define ALPS_UTILITY_SIZE_HPP

#include <alps/type_traits/is_sequence.hpp>

#include <cstddef>
#include <type_traits>

namespace alps {

template <class T>
inline typename std::enable_if<!is_sequence<T>::value,std::size_t>::type
size(T const&) 
{
  return 1;
}

template <class T>
inline typename std::enable_if<is_sequence<T>::value,std::size_t>::type
size(T const& a) 
{
  return a.size();
}

} // end namespace alps

#endif // ALPS_UTILITY_SIZE_HPP
