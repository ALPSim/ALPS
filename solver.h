/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005 by Philipp Werner <werner@comp-phys.org>,
*                       Matthias Troyer <troyer@comp-phys.org>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
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

/* $Id$ */

#ifndef ALPS_DMFT_SOLVER_H
#define ALPS_DMFT_SOLVER_H

/// @file solver.h
/// @brief declares the abstract base class for the impurity solver

#include <alps/parameters.h>
#include <vector>
#include <utility>
#include "types.h"
/// the base class for all impurity solvers
class ImpuritySolver {
public:
  ImpuritySolver() {}
  
  ///  @brief solves the impurity problem
  ///
  ///  @param G0 the inverse Weiss field (without i * omega), as a function of Matsubara 
  ///               frequency omega_n
  ///  @param parms the simulation parameters
  ///  
  ///  Note that  mu, U and h are passed as simulation parameters.
  ///  
  ///  @return the Green's functions for up and down spins as function of Matsubara frequency */
    
  virtual multiple_vector_type solve(const multiple_vector_type& G0, const alps::Parameters& parms)=0;
  virtual ~ImpuritySolver() {}
};

#endif
