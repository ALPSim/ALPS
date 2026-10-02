/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef PARAPACK_FOOTPRINT_H
#define PARAPACK_FOOTPRINT_H

#include <alps/config.h>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>

#include <string>
#include <vector>

namespace alps {

template<typename T>
std::size_t footprint(T const& t, typename boost::enable_if<boost::is_pod<T> >::type* = 0) {
  return sizeof(T);
}

template<typename T>
std::size_t footprint(T const& t, typename boost::disable_if<boost::is_pod<T> >::type* = 0) {
  return t.footprint();
}

template<typename T>
std::size_t footprint(std::vector<T> const& v) {
  return sizeof(std::vector<T>) + sizeof(T) * v.capacity();
}

template<typename C, typename T, typename A>
std::size_t footprint(std::basic_string<C, T, A> const& v) {
  return sizeof(std::basic_string<C, T, A>) + sizeof(C) * v.capacity();
}

} // end namespace alps

#endif // PARAPACK_FOOTPRINT_H
