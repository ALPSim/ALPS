/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Prakash Dayal <prakash@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
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

