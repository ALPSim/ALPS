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
#include <alps/ngs.hpp>
#ifdef ALPS_HAVE_MPI
#include <boost/mpi/environment.hpp>
#endif

int main(int argc, char** argv) {
  alps::mcoptions options(argc, argv);
  if (!options.valid) return 0;
#ifdef ALPS_HAVE_MPI
  boost::mpi::environment env(argc, argv);
#endif
  alps::params parms(alps::hdf5::archive(options.input_file, alps::hdf5::archive::READ));
  try {
    if (options.time_limit != 0)
      throw std::invalid_argument("time limit is passed in the parameter file!");
    if (!parms.defined("MAX_TIME"))
      throw std::runtime_error("parameter MAX_TIME is not defined. How long do you want to run the code for? (in seconds)");
    alps::solvers::cthyb(parms, options.output_file);
  } catch (std::exception const& exc) {
    std::cerr << exc.what() << '\n';
    return -1;
  } catch (...) {
    std::cerr << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
  return 0;
}
