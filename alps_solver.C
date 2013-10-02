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
#include <vector>
#include <utility>
#include <sstream>

#include "boost/assert.hpp"
#include "alps/osiris/comm.h"
#include "alps/scheduler/options.h"
#include "alps/utility/temporary_filename.hpp"

#include "alps_solver.h"
#include "types.h"
#include "boost/filesystem.hpp"
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/exception.hpp>

alps::ImpuritySolver::ImpuritySolver(const scheduler::Factory& factory, int argc, char** argv)
{
	comm_init(argc,argv, true);
	if(is_master()){
          alps::scheduler::NoJobfileOptions opt(1,argv);
          opt.max_check_time=60;
          opt.use_mpi=true;
          master_scheduler = new alps::scheduler::SingleScheduler(opt,factory);
	} else{ //a slave lives for many iterations...
            alps::scheduler::NoJobfileOptions opt(1,argv);
            opt.use_mpi=true;
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
	throw std::logic_error("a slave should never reach this section.");
  }
}


itime_green_function_t  alps::ImpuritySolver::solve(
               const itime_green_function_t & G0,
               alps::Parameters p)
{
  if(!is_master()) throw std::logic_error("a slave should never reach this section.");
  std::string infile = alps::temporary_filename("alps_external_solver_in_");
  std::string outfile = alps::temporary_filename("alps_external_solver_out_");
  if(p.defined("TMPNAME")){
    infile=p["TMPNAME"]+std::string(".in.h5");
    outfile=p["TMPNAME"]+std::string(".out.h5");
  }
  p["INFILE"]=infile;
  p["OUTFILE"]=outfile;

  // write input file
  {
    alps::hdf5::archive solver_input(infile, alps::hdf5::archive::WRITE);
    solver_input<<alps::make_pvp("/parameters",p);
    G0.write_hdf5(solver_input, "/G0");
  }
  int res=solve_it(p);
  if (res)
    boost::throw_exception(
      std::runtime_error(" solver finished with nonzero exit code"));
    
  // now extract the resuls
  dynamic_cast<ImpurityTask*>(get_task())->get_result();
  

  // read the output
  unsigned int N=(unsigned int)p["N"];
  unsigned int sites     =(unsigned int)p.value_or_default("SITES", 1);
  unsigned int flavors   =(unsigned int)p.value_or_default("FLAVORS", 2);
  itime_green_function_t g(N+1, sites, flavors, G0.shape());
  {
    alps::hdf5::archive ar(outfile, alps::hdf5::archive::READ);
    g.read_hdf5(ar, "/G0");
  }
  clear(); // destroy the simulation
  boost::filesystem::remove(boost::filesystem::system_complete(outfile));
  return g;
}

std::pair<matsubara_green_function_t, itime_green_function_t> alps::ImpuritySolver::solve_omega(
                 const matsubara_green_function_t& G0_omega
               , Parameters p)
{
  BOOST_ASSERT(is_master());
  std::string infile = alps::temporary_filename("alps_external_solver_");
  std::string outfile = alps::temporary_filename("alps_external_solver_out_");
  if(p.defined("TMPNAME")){
    infile=p["TMPNAME"]+std::string(".in.h5");
    outfile=p["TMPNAME"]+std::string(".out.h5");
  }
  p["INFILE"]=infile;
  p["OUTFILE"]=outfile;

  {
    alps::hdf5::archive solver_input(infile, alps::hdf5::archive::WRITE);
    solver_input<<alps::make_pvp("/parameters", p);
    G0_omega.write_hdf5(solver_input, "/G0");
  }

  int res=solve_it(p);
  if (res)
    boost::throw_exception(
      std::runtime_error(" solver finished with nonzero exit code"));
    
  // now extract the results
      dynamic_cast<MatsubaraImpurityTask*>(get_task())
        ->get_result();

  // read the output

  unsigned int n_matsubara=(unsigned int)p["NMATSUBARA"];
  unsigned int n_tau=(unsigned int)p["N"];
  unsigned int n_site     =(unsigned int)p.value_or_default("SITES", 1);
  unsigned int n_orbital  =(unsigned int)p.value_or_default("FLAVORS", 2);
  matsubara_green_function_t G_omega(n_matsubara, n_site, n_orbital, G0_omega.shape());
  itime_green_function_t G_tau(n_tau+1, n_site, n_orbital, G0_omega.shape());
  {
    alps::hdf5::archive ar(outfile, alps::hdf5::archive::READ);
    G_omega.read_hdf5(ar, "/G_omega");
    G_tau.read_hdf5(ar, "/G_tau");
  }
 
  clear(); // destroy the simulation
  boost::filesystem::remove(boost::filesystem::system_complete(outfile));
  return std::make_pair(G_omega, G_tau);
}

std::pair<nambu_matsubara_green_function_t, nambu_itime_green_function_t> alps::ImpuritySolver::solve_nambu(
                 const nambu_matsubara_green_function_t& G0_omega
               , Parameters p)
{
  throw std::logic_error("so far we have no NAMBU alps solvers - check if you really want to be here!");
  BOOST_ASSERT(is_master());
  std::string infile = alps::temporary_filename("alps_external_solver_");
  std::string outfile = alps::temporary_filename("alps_external_solver_out_");
  if(p.defined("TMPNAME")){
    infile=p["TMPNAME"]+std::string(".in.h5");
    outfile=p["TMPNAME"]+std::string(".out.h5");
  }
  p["INFILE"]=infile;
  p["OUTFILE"]=outfile;

  {
    alps::hdf5::archive solver_input(infile, alps::hdf5::archive::WRITE);
    solver_input<<alps::make_pvp("/parameters", p);
    G0_omega.write_hdf5(solver_input, "/G0");
  }

  int res=solve_it(p);
  if (res)
    boost::throw_exception(
      std::runtime_error(" solver finished with nonzero exit code"));
    
  // now extract the results
      dynamic_cast<MatsubaraImpurityTask*>(get_task())
        ->get_result();

  // read the output

  unsigned int n_matsubara=(unsigned int)p["NMATSUBARA"];
  unsigned int n_tau=(unsigned int)p["N"];
  unsigned int n_site     =(unsigned int)p.value_or_default("SITES", 1);
  unsigned int n_orbital  =(unsigned int)p.value_or_default("FLAVORS", 2);
  nambu_matsubara_green_function_t G_omega(n_matsubara, n_site, n_orbital, diagonal);
  nambu_itime_green_function_t G_tau(n_tau+1, n_site, n_orbital, diagonal);
  {
    alps::hdf5::archive ar(outfile, alps::hdf5::archive::READ);
    G_omega.read_hdf5(ar, "/G_omega");
    G_tau.read_hdf5(ar, "/G_tau");
  }
 
  clear(); // destroy the simulation
  boost::filesystem::remove(boost::filesystem::system_complete(outfile));
  return std::make_pair(G_omega, G_tau);
}


