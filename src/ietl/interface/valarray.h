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

/* $Id: valarray.h,v 1.4 2003/09/13 10:30:24 troyer Exp $ */

#ifndef IETL_INTERFACE_VALARRAY_H
#define IETL_INTERFACE_VALARRAY_H

#include <ietl/complex.h>
#include <valarray>
#include <numeric>
#include <cmath>

namespace ietl {

  template <class T>
  struct add_abs2 {
    T operator()(T sum, T x) { return sum+x*x;}
  };

  template <class T>
  struct add_abs2<std::complex<T> > {
    T operator()(T sum, const std::complex<T>& x) { return sum+x.real()*x.real()+x.imag()*x.imag();}
  };

  template <class T>
  struct conj_mult {
    T operator()(T x, T y) { return x*y;}
  };

  template <class T>
  struct conj_mult<std::complex<T> > {
    std::complex<T> operator()(const std::complex<T>& x, const std::complex<T>& y) 
    { return std::conj(x)*y; }
  };
  
  
  template <class T>
  T* get_data(std::valarray<T>& c) 
  {
    return &(c[0]);
  }

  template <class T>
  const T* get_data(const std::valarray<T>& c) 
  {
    return get_data(const_cast<std::valarray<T>&>(c));
  }
  
  
  template < class T> 
    void clear(std::valarray<T>& c) {
    std::fill(get_data(c),get_data(c)+c.size(),0.);
  }  

  template < class T, class Gen> 
  void generate(std::valarray<T>& c, const Gen& gen)
  {
    std::generate(get_data(c),get_data(c)+c.size(),gen);
  }

  template < class T> 
  typename real_type<T>::type two_norm(const std::valarray<T>& c)
  {
    return std::sqrt(std::accumulate(get_data(c),get_data(c)+c.size(),0.,add_abs2<T>()));
  }

  template < class T> 
  T dot(const std::valarray<T>& x,const std::valarray<T>& y)
  {
    return std::inner_product(get_data(x),get_data(x)+x.size(),get_data(y),T(),std::plus<T>(),conj_mult<T>());
  }

}

#endif

