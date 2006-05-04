/*******************************************************************************
*
*  DMFT Project
*
* Copyright (C) 2005-2006 by Philipp Werner <werner@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*                            Emanuel Gull <gullc@comp-phys.org>
*
* This software is part of the  Applications, published under the 
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the  Application License along with
* the  Applications; see the file LICENSE.txt. If not, the license is also
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

#include "alps_solver.h"
#include "xml.h"
#include "types.h"
#include <boost/assert.hpp>
#include <alps/osiris/comm.h>
#include <alps/scheduler/options.h>
#include <vector>
#include <utility>
#include <sstream>

alps::ImpuritySolver::ImpuritySolver(const scheduler::Factory& factory, int argc, char** argv)
// : the_scheduler(scheduler::start_single(f,argc,argv))
{
	comm_init(argc,argv);
	if(is_master()){
          alps::scheduler::NoJobfileOptions opt(1,argv);
          master_scheduler = new alps::scheduler::SingleScheduler(opt,factory);
	} else{ //a slave lives for many iterations...
            alps::scheduler::NoJobfileOptions opt(argc,argv);
	    int res=0;
	    alps::scheduler::Scheduler *slave_scheduler = new alps::scheduler::Scheduler(opt,factory);
	    res = slave_scheduler->run();
            delete slave_scheduler;
	    comm_exit();
	    exit(0);
	}
}


alps::ImpuritySolver::~ImpuritySolver()
{
  scheduler::stop_single();
}


int alps::ImpuritySolver::solve_it(Parameters const& p)
{
  if (is_master())
  {
      master_scheduler->create_task(p);
      return master_scheduler->run();
  }
  else{
	std::cerr<<"a slave should never reach this section."<<std::endl;
	exit(1);
  }
  /*std::cout<<"creating task"<<std::endl;
  the_scheduler->create_task(p);
  std::cout<<"running task"<<std::endl;
  return the_scheduler->run();*/
}


itime_green_function_t  alps::ImpuritySolver::solve(
               const itime_green_function_t & G0, 
               const alps::Parameters& p)
{
  BOOST_ASSERT(is_master());
  
  alps::Parameters parms(p);
  
  std::ostringstream G0_text;
  alps::oxstream G0_xml(G0_text);
  write_itime(G0_xml,G0);
  parms["G0"] = G0_text.str();
  
  int res=solve_it(parms);
  if (res)
    boost::throw_exception(
      std::runtime_error(" solver finished with nonzero exit code"));
    
  // now extract the resuls
  itime_green_function_t G = 
      dynamic_cast<ImpurityTask*>(get_task())->get_result();
  
  clear(); // destroy the simulation
  return G;
}

std::pair<matsubara_green_function_t, itime_green_function_t>
alps::ImpuritySolver::solve_omega(
                 const matsubara_green_function_t& G0_omega
               , const Parameters& p)
{
  BOOST_ASSERT(is_master());
  
  alps::Parameters parms(p);
  std::ostringstream G0_omega_text;
  alps::oxstream G0_omega_xml(G0_omega_text);
  
  write_freq(G0_omega_xml,G0_omega);
  parms["G0(omega)"] = G0_omega_text.str();
  int res=solve_it(parms);
  if (res)
    boost::throw_exception(
      std::runtime_error(" solver finished with nonzero exit code"));
    
  // now extract the results
  std::pair<matsubara_green_function_t, itime_green_function_t> G = 
      dynamic_cast<MatsubaraImpurityTask*>(get_task())
        ->get_result();

  clear(); // destroy the simulation
  return G;

}



