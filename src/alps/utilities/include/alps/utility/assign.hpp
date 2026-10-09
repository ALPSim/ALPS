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

#ifndef ALPS_UTILITY_ASSIGN_HPP
#define ALPS_UTILITY_ASSIGN_HPP

#include <valarray>

namespace alps {

template <class X, class Y> 
inline void assign(X& x,const Y& y) 
{
  x=y;
}

template <class X, class Y> 
inline void assign(std::valarray<X>& x, std::valarray<Y> const& y) 
{
  x.resize(y.size()); 
  for (std::size_t i=0;i<y.size();++i) 
    x[i]=y[i];
}

template <class X> 
inline void assign(std::valarray<X>& x, std::valarray<X> const& y) 
{
  x.resize(y.size()); 
  x=y;
}

} // end namespace alps

#endif // ALPS_UTILITY_ASSIGN_HPP
