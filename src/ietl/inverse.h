/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2002 by Rene Villiger <rvilliger@smile.ch>,
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

/* $Id: inverse.h,v 1.9 2004/02/15 23:30:42 troyer Exp $ */

#ifndef IETL_INVERSE_H
#define IETL_INVERSE_H

#include <ietl/traits.h>
#include <ietl/complex.h>

// REQUIRES class SOLVER with interface
//
//    void operator()(MATRIX matrix, MAGNITUDE_TYPE sigma, VECTOR_TYPE v, VECTOR_TYPE y)
// 
// Solves the equation
//    y = (A - sigma*I)^{-1} * v
//

namespace ietl
{
   template <class MATRIX, class GEN, class SOLVER, class ITER, class VS>
   std::pair<typename vectorspace_traits<VS>::magnitude_type,
             typename vectorspace_traits<VS>::vector_type>
      inverse(const MATRIX& matrix, GEN& gen,const SOLVER& solver, ITER& iter,
                          typename vectorspace_traits<VS>::magnitude_type sigma,
              const VS& vec)
      {
         typedef typename vectorspace_traits<VS>::vector_type vector_type;
         typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
         typedef typename vectorspace_traits<VS>::magnitude_type magnitude_type;
         
         vector_type y   = new_vector(vec);
         vector_type v   = new_vector(vec);
         magnitude_type theta;
         magnitude_type residual;

         // Start with vector y=z, the initial guess
         ietl::generate(y,gen);
         ietl::project(y,vec);
         
         // start iteration loop
         do {
            v=(1./ietl::two_norm(y)) * y;
            try {
                          solver(matrix, sigma, v, y);  // y = (A-\sigma*I)^{-1} v
            }
            catch (...) {
              break; // done with iteration
            }  
                        
            theta = ietl::real(ietl::dot(v,y));
            v = y-theta*v;  // residual = | y - \theta*v |_2
            residual = ietl::two_norm(v);
            ++iter;
            
         // check for convergence
         } while(!iter.finished(residual, theta));

         // accept \lambda = \sigma + 1/\theta  and  x = y/\theta 
         y/=theta;
         theta = sigma + 1./theta;
         
         return std::make_pair(theta, y);
      }
}

#endif
