/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Prakash Dayal <prakash@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: blitz.h,v 1.10 2003/09/05 08:12:38 troyer Exp $ */

#ifndef IETL_INTERFACE_BLITZ_H
#define IETL_INTERFACE_BLITZ_H

#include <ietl/complex.h>
#include <blitz/array.h>
#include <blitz/vecglobs.h>
#include <ietl/traits.h>


namespace ietl {
  
  template < class Cont>
    void clear(Cont& c) {
    std::fill(c.begin(),c.end(),0.);
  }

  template < class Cont, class Gen> 
    void generate(Cont& c, const Gen& gen) {
    std::generate(c.begin(),c.end(),gen);
  }

  template <class T, int D>
  typename number_traits<T>::magnitude_type two_norm(const blitz::Array<T,D>& v) {
    return std::sqrt(ietl::real(dot(v,v)));
  }

  template <class T, int D>
  T dot(const blitz::Array<T,D>& x, const blitz::Array<T,D>& y) {
    return blitz::sum(x*y);
  }

  template <class T, int D>
  T dot(const blitz::Array<std::complex<T>,D>& x, const blitz::Array<std::complex<T>,D>& y) {
    
    return blitz::sum(blitz::conj(x)*y);
  }

  template <class T, int D>
  void copy(const blitz::Array<T,D>& x, blitz::Array<T,D>& y) {
    y=x;
    y.makeUnique();
  }
}

namespace std {  
  template <class T, int D>
  void swap(blitz::Array<T,D>& x, blitz::Array<T,D>& y) {
    blitz::cycleArrays(x,y);
  }

}

#endif // IETL_INTERFACE_BLITZ_H

