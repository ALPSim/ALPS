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

/* $Id: alps_solver.C 342 2009-01-28 22:31:54Z fuchs $ */

#include "alps_solver.h"
#include "xml.h"
#include "types.h"
#include <boost/assert.hpp>
#include <alps/osiris/comm.h>
#include <alps/scheduler/options.h>
#include <vector>
#include <utility>
#include <sstream>

#define DEFAULT_CHECK_TIME 300

alps::ImpuritySolver::ImpuritySolver(const scheduler::Factory& factory, int argc, char** argv, bool h5input)
{
#ifdef ALPS_HAVE_MPI
	comm_init(argc,argv, true);
#else
	comm_init(argc,argv, false);
#endif
	if(is_master()){
    alps::scheduler::NoJobfileOptions opt(1,argv);
    unsigned int max_time;
    if (!h5input){
      alps::Parameters parms;
      std::ifstream is(argv[1]);
      is>>parms;
      max_time = (unsigned int)(parms.value_or_default("MAX_TIME",DEFAULT_CHECK_TIME));
    } else {
      alps::Parameters parms;
      alps::hdf5::archive ar(argv[1], "r");
      ar["/parameters"] >> parms;
      max_time = (unsigned int)(parms.value_or_default("MAX_TIME",DEFAULT_CHECK_TIME));
    }
    if (max_time < DEFAULT_CHECK_TIME) { opt.checkpoint_time=max_time; opt.max_check_time=max_time; opt.min_check_time=max_time; }
    else {
      unsigned int checks = (max_time + DEFAULT_CHECK_TIME - 1) / DEFAULT_CHECK_TIME;
      opt.checkpoint_time=(max_time+checks-1)/checks;
      opt.max_check_time=opt.checkpoint_time; opt.min_check_time=opt.checkpoint_time; 
    }
    master_scheduler = new alps::scheduler::SingleScheduler(opt,factory);
	} else{ //a slave lives for many iterations...
    alps::scheduler::NoJobfileOptions opt(1,argv);
    alps::scheduler::Scheduler *slave_scheduler = new alps::scheduler::Scheduler(opt,factory);
    slave_scheduler->run();
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
    abort();
  }
}


itime_green_function_t  alps::ImpuritySolver::solve(const itime_green_function_t & G0, const alps::Parameters& p)
{
  BOOST_ASSERT(is_master());
  
  alps::Parameters parms(p);
  std::string basename=parms["BASENAME"];
  //boost::filesystem::remove(basename+".h5");
  alps::hdf5::archive dumpfile(basename+".h5", "a");
  dumpfile["/parameters"]<<parms;
  
  std::ostringstream G0_text;
  alps::oxstream G0_xml(G0_text);
  write_itime(G0_xml,G0);
  parms["G0"] = G0_text.str();
  
  int res=solve_it(parms);
  //boost::filesystem::remove(basename+".h5");
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
alps::ImpuritySolver::solve_omega(const matsubara_green_function_t& G0_omega, const Parameters& p)
{
  BOOST_ASSERT(is_master());
  alps::Parameters parms(p);
  
  std::string basename=parms["BASENAME"];
  //boost::filesystem::remove(basename+".h5");
  alps::hdf5::archive dumpfile(basename+".h5", "a");
  dumpfile["/parameters"]<<parms;
  
  std::ostringstream G0_omega_text;
  alps::oxstream G0_omega_xml(G0_omega_text);
  
  write_freq(G0_omega_xml,G0_omega);
  parms["G0(omega)"] = G0_omega_text.str();
  int res=solve_it(parms);
  //boost::filesystem::remove(basename+".h5");
  if (res)
    boost::throw_exception(
                           std::runtime_error(" solver finished with nonzero exit code"));
  
  std::pair<matsubara_green_function_t, itime_green_function_t> G = 
  dynamic_cast<MatsubaraImpurityTask*>(get_task())
  ->get_result();
  clear(); // destroy the simulation
  return G;
  
}



