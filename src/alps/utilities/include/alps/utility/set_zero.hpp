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

#ifndef ALPS_UTILITY_SET_ZERO_HPP
#define ALPS_UTILITY_SET_ZERO_HPP

#include <alps/type_traits/is_sequence.hpp>
#include <alps/type_traits/element_type.hpp>

#include <boost/utility/enable_if.hpp>

#include <algorithm>

namespace alps {

template <class X> 
inline typename boost::disable_if<is_sequence<X>,void>::type
set_zero(X& x) { x=X();}

template <class X> 
inline typename boost::enable_if<is_sequence<X>,void>::type
set_zero(X& x) 
{
  std::fill(x.begin(),x.end(),typename element_type<X>::type());
}

//Enable set_zero for valarrays
template <class X> 
inline typename boost::enable_if<is_sequence<std::valarray<X> >,void>::type
set_zero(std::valarray<X> & x) 
{
  for(int i=0;i<x.size();++i) x[i]=X();
}


} // end namespace alps

#endif // ALPS_UTILITY_RESIZE_HPP
