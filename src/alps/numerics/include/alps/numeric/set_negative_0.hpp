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

#ifndef ALPS_NUEMRIC_SET_NEGATIVE_0_HPP
#define ALPS_NUEMRIC_SET_NEGATIVE_0_HPP

#include <alps/type_traits/is_sequence.hpp>
#include <boost/utility/enable_if.hpp>
#include <complex>


namespace alps { namespace numeric {

template <class T>
inline typename boost::disable_if<is_sequence<T>,void>::type
set_negative_0(T& x)
{
  if (x<T()) 
    x=T();
}

template <class T>
inline void set_negative_0(std::complex<T>& x)
{ 
  if (std::real(x)<0. || std::imag(x)<0.) 
    x=0.;
}

template <class T>
inline typename boost::enable_if<is_sequence<T>,void>::type
set_negative_0(T& a) 
{
  for(std::size_t i=0; i!=a.size(); ++i)
    set_negative_0(a[i]);
}



} } // end namespace alps::numeric

#endif // ALPS_NUEMRIC_SET_NEGATIVE_0_HPP
