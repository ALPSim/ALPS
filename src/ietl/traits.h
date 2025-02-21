/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2002 by Prakash Dayal <prakash@comp-phys.org>,
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

/* $Id: traits.h,v 1.2 2004/02/15 23:30:42 troyer Exp $ */

#ifndef IETL_VECTORSPACETRAITS__H
#define IETL_VECTORSPACETRAITS__H

#include <complex>

namespace ietl {
  template <class T> 
    struct number_traits {
      typedef T magnitude_type;
    };
    
  template <class T>
    struct number_traits<std::complex<T> > {
      typedef T magnitude_type;
    };

  template <class VS>
    struct vectorspace_traits {
      typedef typename VS::vector_type vector_type;
      typedef typename VS::size_type size_type;
      typedef typename VS::scalar_type scalar_type;
          typedef typename number_traits<scalar_type>::magnitude_type magnitude_type;
    };
    
}
#endif
