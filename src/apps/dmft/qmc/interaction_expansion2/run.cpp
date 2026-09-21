/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2010 by Emanuel Gull <gull@phys.columbia.edu>,
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

#include <alps/solvers.hpp>
#include "interaction_expansion.hpp"
#include "fouriertransform.h"
#include <alps/utility/copyright.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#ifdef ALPS_HAVE_MPI
#include <alps/mcmpiadapter.hpp>
using sim_type = alps::mcmpiadapter<HubbardInteractionExpansion>;
#else
using sim_type = HubbardInteractionExpansion;
#endif


void compute_greens_functions(const alps::results_type<HubbardInteractionExpansion>::type&,
                             alps::params const&, std::string const&);

namespace {
bool stop_requested(boost::posix_time::ptime const & end_time) {
//stops the simulation if time > end_time or if signals received.
  static alps::ngs::signal signal;
  return !signal.empty() || boost::posix_time::second_clock::local_time() > end_time;
}
}

void alps::solvers::ctint(alps::params const& parms, std::string const& output_file) {
  int rank;
#ifndef ALPS_HAVE_MPI
  rank=0;
  sim_type s(parms,rank);
#else
  boost::mpi::communicator c;
  c.barrier();
  rank=c.rank();
  sim_type s(parms, c);
#endif
  if (rank==0) {
    alps::print_copyright(std::cout);
    std::cout << "****************************************************************"<<std::endl;
    std::cout << "* Please cite in scientific publications:                      *"<<std::endl;
    std::cout << "* We used the ALPS [1] implementation [2] of the CT-INT        *"<<std::endl;
    std::cout << "* interaction expansion CT-QMC solver [3,4,5].                 *"<<std::endl;
    std::cout << "* [1] JSTAT (2011) P05001; [2] CPC 182, 1078 (2011);           *"<<std::endl;
    std::cout << "* [3] JETP Lett. 80, 61 (2004); [4] PRB 72, 035122 (2005);     *"<<std::endl;
    std::cout << "* [5] RMP 83, 349 (2011).                                      *"<<std::endl;
    std::cout << "****************************************************************"<<std::endl;
  }
  //run the simulation
  s.run(boost::bind(&stop_requested, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds((int)parms["MAX_TIME"])));

  //on the master: collect MC results and store them in file, then postprocess
  if (rank==0){
    alps::results_type<HubbardInteractionExpansion>::type results = collect_results(s);
    save_results(results, parms, output_file, "/simulation/results");
    //compute the output Green's function and Fourier transform it, store in the right path
    compute_greens_functions(results, parms, output_file);
  }
#ifdef ALPS_HAVE_MPI
  else{
  collect_results(s);
  }
  c.barrier();
#endif
}
