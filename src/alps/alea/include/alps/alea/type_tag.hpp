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

#ifndef ALPS_ALEA_TYPE_TAG_H
#define ALPS_ALEA_TYPE_TAG_H

#include <alps/type_traits/type_tag.hpp>
#include <vector>
#include <valarray>

namespace alps {

template <class T>
struct type_tag<std::valarray<T> >
 : public boost::mpl::int_<256 + type_tag<T>::value> {};

template <class T>
struct type_tag<std::vector<T> >
 : public boost::mpl::int_<256 + type_tag<T>::value> {};

} // end namespace alps

#endif // ALPS_ALEA_TYPE_TAG_H
