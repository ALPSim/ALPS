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
#include <alps/ngs/mcoptions.hpp>
#include <alps/params.hpp>
#include <alps/parameter_input.hpp>
#include <iostream>
#include <stdexcept>
#ifdef ALPS_HAVE_MPI
#include <alps/utilities/mpi.hpp>
#endif

int main(int argc, char** argv) {
#ifdef ALPS_HAVE_MPI
  // outside the try: destroyed while an exception unwinds, it calls MPI_Abort
  // before the catch below can print the message. The catch aborts itself,
  // since other ranks may be waiting in a collective.
  alps::mpi::environment env(argc, argv);
#endif
  try {
    alps::mcoptions options(argc, argv, "interaction");
    if (!options.valid) return 0;
    alps::params parms = alps::read_params_file(options.input_file);
    if (options.time_limit != 0)
      throw std::invalid_argument("time limit is passed in the parameter file!");
    if (!parms.supplied("MAX_TIME"))
      throw std::runtime_error("parameter MAX_TIME is not defined. How long do you want to run the code for? (in seconds)");
    alps::solvers::ctint(parms, options.output_file);
  } catch (std::exception const& exc) {
    std::cerr << exc.what() << '\n';
#ifdef ALPS_HAVE_MPI
    alps::mpi::environment::abort(255);
#endif
    return -1;
  } catch (...) {
    std::cerr << "Fatal Error: Unknown Exception!\n";
#ifdef ALPS_HAVE_MPI
    alps::mpi::environment::abort(255);
#endif
    return -2;
  }
  return 0;
}
