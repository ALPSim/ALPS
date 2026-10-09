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

#ifndef ALPS_NUMERIC_UPDATE_MINMAX_HPP
#define ALPS_NUMERIC_UPDATE_MINMAX_HPP

#include <alps/type_traits/slice.hpp>

namespace alps { namespace numeric {

template <class T>
void update_max(T& lhs, T const& rhs)
{
  for (typename slice_index<T>::type it = slices(lhs).first; 
       it < slices(lhs).second && it < slices(rhs).second; ++it)
    if (slice_value(lhs,it) < slice_value(rhs,it))
      slice_value(lhs,it) = slice_value(rhs,it);
}

template <class T>
void update_min(T& lhs, T const& rhs)
{
  for (typename slice_index<T>::type it = slices(lhs).first; 
       it < slices(lhs).second && it < slices(rhs).second; ++it)
    if (slice_value(rhs,it) < slice_value(lhs,it))
      slice_value(lhs,it) = slice_value(rhs,it);
}



} } // end namespace alps::numeric

#endif // ALPS_NUMERIC_UPDATE_MINMAX_HPP
