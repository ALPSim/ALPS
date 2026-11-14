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

// Usage: ct-hyb-segment parameters.ini [NAME=value ...]
//        ct-hyb-segment --help
// Results are written to sim.h5 in the working directory.
#include <alps/solvers.hpp>
#include <alps/params.hpp>
#include "hyb.hpp"
#include <iostream>
#include <stdexcept>
#ifdef ALPS_HAVE_MPI
#include <alps/utilities/mpi.hpp>
#endif

int main(int argc, char* argv[]){
#ifdef ALPS_HAVE_MPI
  // outside the try: destroyed while an exception unwinds, it calls MPI_Abort
  // before the catch below can print the message. The catch aborts itself,
  // since other ranks may be waiting in a collective.
  alps::mpi::environment env(argc, argv);
#endif
  try {
    alps::params parameters(argc, (const char**)argv, "/parameters");
    cthyb_segment::hybridization::define_parameters(parameters);
    if (parameters.help_requested(std::cout))
      return 0;
    alps::solvers::cthyb_segment(parameters, "sim.h5");
  }
  catch(std::exception& exc){
    std::cerr<<exc.what()<<std::endl;
#ifdef ALPS_HAVE_MPI
    alps::mpi::environment::abort(255);
#endif
    return -1;
  }
  catch(...){
    std::cerr << "Fatal Error: Unknown Exception!\n";
#ifdef ALPS_HAVE_MPI
    alps::mpi::environment::abort(255);
#endif
    return -2;
  }
  return 0;
}
