/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2002 by Prakash Dayal <prakash@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: complex.h,v 1.5 2003/09/13 10:30:24 troyer Exp $ */

#ifndef IETL_COMPLEX_H
#define IETL_COMPLEX_H 

#include <complex>

namespace ietl {
  template <class T>
  struct real_type {
    typedef T type;
  };
  
  template <class T>
  struct real_type<std::complex<T> > {
    typedef T type;
  };
  
  template <class T> T real (T x) { return x;}
  template <class T> T real (std::complex<T> x) { return x.real();}
  template <class T> T conj (T x) { return x;}
  template <class T> T conj (std::complex<T> x) { return std::conj(x);}
}
#endif

