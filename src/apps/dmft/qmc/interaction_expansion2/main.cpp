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
#include <alps/hdf5/archive.hpp>
#include <alps/ngs/mcoptions.hpp>
#include <alps/params.hpp>
#include <iostream>
#include <stdexcept>
#ifdef ALPS_HAVE_MPI
#include <alps/utilities/mpi.hpp>
#endif

int main(int argc, char** argv) {
  try {
    alps::mcoptions options(argc, argv, "interaction");
    if (!options.valid) return 0;
#ifdef ALPS_HAVE_MPI
    alps::mpi::environment env(argc, argv);
#endif
    alps::params parms;
    {
      alps::hdf5::archive input(options.input_file, alps::hdf5::archive::READ);
      input["/parameters"] >> parms;
    }
    if (options.time_limit != 0)
      throw std::invalid_argument("time limit is passed in the parameter file!");
    if (!parms.supplied("MAX_TIME"))
      throw std::runtime_error("parameter MAX_TIME is not defined. How long do you want to run the code for? (in seconds)");
    alps::solvers::ctint(parms, options.output_file);
  } catch (std::exception const& exc) {
    std::cerr << exc.what() << '\n';
    return -1;
  } catch (...) {
    std::cerr << "Fatal Error: Unknown Exception!\n";
    return -2;
  }
  return 0;
}
