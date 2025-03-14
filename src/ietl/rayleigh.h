/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Rene Villiger <rvilliger@smile.ch>,
*                            Prakash Dayal <prakash@comp-phys.org>,
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

/* $Id: rayleigh.h,v 1.13 2004/02/15 23:30:42 troyer Exp $ */

#ifndef IETL_RAYLEIGH_H
#define IETL_RAYLEIGH_H

#include <ietl/traits.h>
#include <ietl/complex.h>

// REQUIRES class SOLVER
//
//    bool solver(MATRIX matrix, MAGNITUDE_TYPE rho, VECTOR_TYPE v, VECTOR_TYPE y)
// 
// Solves the equation
//    y = (A - rho*I)^{-1} * v
//
// returns true, if singular, false otherwise.
// 

namespace ietl
{
  template <class MATRIX, class GEN, class SOLVER, class ITER, class VS>
  std::pair<typename ietl::number_traits<typename vectorspace_traits<VS>::scalar_type>::magnitude_type,
            typename vectorspace_traits<VS>::vector_type>
  rayleigh(const MATRIX& matrix, 
               GEN& gen, 
               SOLVER& solver, 
               ITER& iter,
               const VS& vec )
  {
    typedef typename vectorspace_traits<VS>::vector_type vector_type;
    typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
    typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
         
    vector_type y = new_vector(vec);
    vector_type v = new_vector(vec);
    magnitude_type rho;
    magnitude_type theta;

    // v = y / |y|_2 and \rho_1 = \rho(v)
    ietl::generate(y,gen);
    ietl::project(y,vec);
    v = (1./ietl::two_norm(y))*y;
    ietl::mult(matrix, v, y);
    rho = ietl::real(ietl::dot(y, v) / ietl::two_norm(v));
    
    // start iteration
    do {
      // y = (A - \rho_k I)^{-1} v  (if singular stop iteration)
      try {
        solver(matrix, rho, v, y);
      }
      catch (...) {
        break; // done with iteration
      }
        
      theta = ietl::two_norm(y);
            
      // \rho_{k+1} = \rho_k + y^\star * v / theta^2
      rho += ietl::real(ietl::dot(y,v) / (theta*theta));
      v=(1./theta)*y;
      ++iter;
      // if \theta > \varepsilon_M ^{-1/2}, stop
    }  while(!iter.finished(1./(theta*theta),0.));
         
    // accept \lambda = \rho_k for most recent k and x=v
    return std::make_pair(rho, v);
  }
}

#endif
            
         
