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
