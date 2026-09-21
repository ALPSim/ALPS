/****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2012 by Emanuel Gull <gull@pks.mpg.de>,
 *                       Hartmut Hafermann <hafermann@cpht.polytechnique.fr>
 *
 *  based on an earlier version by Philipp Werner and Emanuel Gull
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include <alps/solvers.hpp>
#include "hyb.hpp"
#include "hybevaluate.hpp"
#include <alps/utility/copyright.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#ifdef ALPS_HAVE_MPI
#include <alps/mcmpiadapter.hpp>
using sim_type = alps::mcmpiadapter<hybridization>;
#else
using sim_type = hybridization;
#endif


int global_mpi_rank;

namespace {
bool stop_requested(boost::posix_time::ptime const & end_time) {
  static alps::ngs::signal signal;
  return !signal.empty() || boost::posix_time::second_clock::local_time() > end_time;
}
void master_final_tasks(const alps::results_type<hybridization>::type &results,
                        const alps::parameters_type<hybridization>::type &parms,
                        const std::string &output_name){
  //do some post processing: collect Green functions and write
  //them into hdf5 files; calls compute vertex at the very end

  alps::hdf5::archive solver_output(output_name, "a");

  evaluate_basics(results,parms,solver_output);
  evaluate_time(results,parms,solver_output);
  evaluate_freq(results,parms,solver_output);
  evaluate_legendre(results,parms,solver_output);
  evaluate_nnt(results,parms,solver_output);
  evaluate_nnw(results,parms,solver_output);
  evaluate_sector_statistics(results,parms,solver_output);
  evaluate_2p(results, parms, solver_output);
}

}

void alps::solvers::cthyb(alps::params const& parms, std::string const& output_file) {
#ifndef ALPS_HAVE_MPI
  global_mpi_rank=0;
  sim_type s(parms,global_mpi_rank);
#else
  boost::mpi::communicator c;
  c.barrier();
  global_mpi_rank=c.rank();
  sim_type s(parms, c);
#endif
  if (global_mpi_rank==0) {
    alps::print_copyright(std::cout);
    std::cout << "****************************************************************"<<std::endl;
    std::cout << "* Please cite in scientific publications:                      *"<<std::endl;
    std::cout << "* We used the ALPS [1] implementation [2] of the `segment'     *"<<std::endl;
    std::cout << "* CT-HYB solver [3,4,5].                                       *"<<std::endl;
    std::cout << "* [1] JSTAT (2011) P05001; [2] CPC 182, 1078 (2011);           *"<<std::endl;
    std::cout << "* [3] PRL 97, 076405 (2006); [4] RMP 83, 349 (2011);           *"<<std::endl;
    std::cout << "* [5] PRB 84, 075145 (2011).                                   *"<<std::endl;
    std::cout << "****************************************************************"<<std::endl;
  }
  //run the simulation
  s.run(boost::bind(&stop_requested, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds((int)parms["MAX_TIME"])));

  //on the master: collect MC results and store them in file, then postprocess
  if (global_mpi_rank==0){
    alps::results_type<hybridization>::type results = collect_results(s);
    std::string output_path = boost::lexical_cast<std::string>(parms["BASEPATH"]|"")+"/simulation/results";
    save_results(results, parms, output_file, output_path); //"/simulation/results");
    master_final_tasks(results, parms, output_file);
  }
#ifdef ALPS_HAVE_MPI
  else{ //on any slave: send back results to master.
    collect_results(s);
  }
  c.barrier();
#endif
}
