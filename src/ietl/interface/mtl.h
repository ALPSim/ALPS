/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2004 by Prakash Dayal <prakash@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: mtl.h,v 1.10 2004/02/11 08:10:09 troyer Exp $ */

#ifndef IETL_MTL_H
#define IETL_MTL_H

#error The MTL-2 interface is not yet implemented

namespace ietl {
  // generate function.     
  template < class Cont, class Gen> 
    void generate(Cont& c, const Gen& gen) {
    std::generate(c.begin(),c.end(),gen);
  }  

  template < class Cont> 
    void clear(Cont& c) {
    std::fill(c.begin(),c.end(),0.);
  }  

}
#endif
