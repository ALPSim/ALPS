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

/* $Id: power.h,v 1.12 2004/02/15 23:30:42 troyer Exp $ */

#ifndef IETL_POWER_H
#define IETL_POWER_H

#include <ietl/traits.h>

namespace ietl {
  template <class MATRIX, class GEN, class IT, class VS>
    std::pair<typename vectorspace_traits<VS>::scalar_type, typename vectorspace_traits<VS>::vector_type>
    power(const MATRIX& m, GEN& start, IT& iter, const VS& vec) {      

    typedef typename vectorspace_traits<VS>::vector_type vector_type;
    typedef typename vectorspace_traits<VS>::scalar_type scalar_type; 
    typedef typename vectorspace_traits<VS>::magnitude_type magnitude_type; 
    
    vector_type vec1 = new_vector(vec);
    ietl::generate(vec1,start);
    project(vec1,vec);
    vector_type vec2 = new_vector(vec);     
    scalar_type lambda;
    magnitude_type residual; 

    do {
      ietl::mult(m,vec1,vec2);
      lambda = ietl::dot(vec1,vec2);
      vec1 *= -lambda;
      vec1+=vec2;
      residual = ietl::two_norm(vec1); 
      vec1=(1./ietl::two_norm(vec2))*vec2;  
      ++iter;
    } while(!iter.finished(residual,lambda));    
    return std::make_pair(lambda,vec1);
  }         
}

#endif
