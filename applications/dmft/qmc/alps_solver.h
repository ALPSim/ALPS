 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>
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

/* $Id: alps_solver.h 338 2009-01-20 01:22:05Z fuchs $ */

#ifndef ALPS_DMFT_ALPS_SOLVER_H
#define ALPS_DMFT_ALPS_SOLVER_H

/// @file solver.h
/// @brief defines the ALPS solvers

#include <alps/scheduler.h>
#include <vector>
#include <utility>
#include "types.h"
#include "solver.h"

namespace alps {

/// the task base class for parallel ALPS solvers

/// the base class for (potentially parallel) solvers using the alps single 
/// scheduler
///
/// Please note that only one such solver may exist. If you have more solvers, 
/// please use a combined factory.

class ImpuritySolver : public ::ImpuritySolver, public ::MatsubaraImpuritySolver
{
public:
  ImpuritySolver(const scheduler::Factory& f, int argc=0, char** argv=0, bool h5input=false);
  ~ImpuritySolver();
    
  scheduler::AbstractTask* get_task() const 
  { 
    return master_scheduler ? master_scheduler->get_task() : 0;
  }
  void checkpoint(const boost::filesystem::path &/*fn*/){
    if(master_scheduler) master_scheduler->checkpoint();
  }
  void clear()
  {
    if (master_scheduler)
      master_scheduler->destroy_task();
  }
  
  itime_green_function_t solve(
    const itime_green_function_t& G0, 
    const Parameters& parms =Parameters());
    
  std::pair<matsubara_green_function_t, itime_green_function_t> solve_omega(
      const matsubara_green_function_t& G0_omega
    , const Parameters& parms=Parameters());

protected:
  int solve_it(Parameters const& p);
  scheduler::SingleScheduler* master_scheduler;
  //scheduler::SingleScheduler* master_scheduler;
  int argc_;
  char **argv_;
};


class ImpurityTask
{
public:
  ImpurityTask() {}
  virtual ~ImpurityTask() {}
  virtual itime_green_function_t get_result() const=0;
};

class MatsubaraImpurityTask
{
public:
  MatsubaraImpurityTask() {}
  virtual ~MatsubaraImpurityTask() {}
  virtual std::pair<matsubara_green_function_t,itime_green_function_t> get_result() =0;
};


}

#endif
