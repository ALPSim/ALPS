 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
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

/* $Id: solver.h 242 2007-03-22 19:44:45Z fuchs $ */

#ifndef ALPS_DMFT_SOLVER_H
#define ALPS_DMFT_SOLVER_H

/// @file solver.h
/// @brief declares the abstract base class for the impurity solver

#include "types.h"
#include "green_function.h"

#include <alps/parameter.h>
#include <vector>
#include <utility>

/// the base class for all impurity solvers taking a Green's function in imaginary time and returning one also in imaginary time.
class ImpuritySolver {
public:
  ImpuritySolver() {}
  
  typedef itime_green_function_t result_type;
  
  ///  @brief solves the impurity problem
  ///
  ///  @param G0 the bare Green's function in imaginary time, 
  ///            provided as an input
  ///              
  ///  @param parms the simulation parameters
  ///               Note that  mu, U and h are passed as simulation parameters 
  ///               inside parms.
  ///  
  ///  @return the Green's functions for up and down spins as function of 
  ///          imaginary time 
    
  virtual result_type solve(
              const itime_green_function_t& G0
            , const alps::Parameters& parms)=0;
            
  virtual ~ImpuritySolver() {}
};

/// the base class for all impurity solvers taking a Green's function in imaginary 
/// time and Matsubara frequency and returning one in Matsubara frequency.
class MatsubaraImpuritySolver{
public:
  MatsubaraImpuritySolver() {}
    
  typedef std::pair<matsubara_green_function_t,itime_green_function_t> result_type;
    
    ///  @brief solves the impurity problem taking both G0_tau and G0_omega
    ///
    ///  @param G0_omega the bare Green's function in Matsubara frequency
    ///  @param G0_tau the bare Green's function in imaginary time
    ///  @param parms the simulation parameters
    ///  
    ///  Note that  mu, U and h are passed as simulation parameters inside parms.
    /// 
    ///  @return the Green's functions for up and down spins as function of 
    ///          Matsubara frequency
    
  virtual result_type solve_omega(const matsubara_green_function_t& G0_omega,
                  const alps::Parameters& parms =alps::Parameters())=0;
            
  virtual ~MatsubaraImpuritySolver() {}
};



#endif
