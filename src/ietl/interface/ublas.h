/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Prakash Dayal <prakash@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/* $Id: ublas.h,v 1.4 2003/09/05 09:27:53 prakash Exp $ */

#ifndef IETL_UBLAS_H
#define IETL_UBLAS_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <ietl/traits.h>

namespace ietl {

  template < class T, class Gen> 
    inline void generate(boost::numeric::ublas::vector<T>& c, Gen& gen) {
    std::generate(c.begin(),c.end(),gen);
  }  

  template < class T> 
    inline void clear(boost::numeric::ublas::vector<T>& c) {
    c.clear();
  }  


  template < class T, class S>
  inline T dot(const boost::numeric::ublas::vector<T,S>& x , const boost::numeric::ublas::vector<T,S>& y) {
   return boost::numeric::ublas::inner_prod (boost::numeric::ublas::conj(x), y);
  }

  template < class T>
   inline typename number_traits<T>::magnitude_type two_norm(boost::numeric::ublas::vector<T>& x) {
   return boost::numeric::ublas::norm_2(x);
  }

  template < class T>
  void copy(const boost::numeric::ublas::vector<T>& x,boost::numeric::ublas::vector<T>& y) {
    y.assign(x);
  }

  template <class M, class T>
  inline void mult(M& m, const boost::numeric::ublas::vector<T>& x, boost::numeric::ublas::vector<T>& y) {
   y=boost::numeric::ublas::prod(m,x);
 }
}

#endif
