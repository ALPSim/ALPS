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

/* $Id: blas.h,v 1.10 2004/02/11 08:10:09 troyer Exp $ */

#error The BLAS interface is not yet implemented 
namespace ietl {
  template < class Cont, class Gen> 
  void generate(Cont& c, const Gen& gen) {
    std::generate(itl::get_data(c),itl::get_data(c)+c.size(),gen);
  }

}
