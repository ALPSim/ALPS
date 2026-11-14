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
//          results are written to sim.h5 in the working directory
//        ct-hyb-segment input.h5 output.h5
//          parameters (and input data) from an HDF5 or ALPS parameter file
//          written by the ALPS DMFT loop; results go to output.h5
//        ct-hyb-segment --help
#include <alps/solvers.hpp>
#include <alps/params.hpp>
#include <alps/parameter_input.hpp>
#include <alps/cli.hpp>
#include "hyb.hpp"
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#ifdef ALPS_HAVE_MPI
#include <alps/utilities/mpi.hpp>
#endif

namespace {
bool is_hdf5_file(const char* name) {
  static const char signature[] = {char(137), 'H', 'D', 'F', '\r', '\n', char(26), '\n'};
  char first[sizeof(signature)];
  std::ifstream file(name, std::ios::binary);
  return file.read(first, sizeof(first)) && std::memcmp(first, signature, sizeof(first)) == 0;
}

// The ALPS external-solver convention: exactly an input file and an output file.
bool is_input_output_call(int argc, char* argv[]) {
  return argc == 3 && argv[1][0] != '-' && !std::strchr(argv[1], '=')
      && argv[2][0] != '-' && !std::strchr(argv[2], '=') && is_hdf5_file(argv[1]);
}
}

int main(int argc, char* argv[]){
  if (alps::handle_cli_information(argc, argv, "cthyb_segment", [&] {
        std::cout << "Usage: " << argv[0] << " parameters.ini [NAME=value ...]\n"
                  << "       " << argv[0] << " input.h5 output.h5\n\n";
        alps::params parameters;
        cthyb_segment::hybridization::define_parameters(parameters);
        parameters.print_help(std::cout);
      }))
    return 0;
#ifdef ALPS_HAVE_MPI
  // outside the try: destroyed while an exception unwinds, it calls MPI_Abort
  // before the catch below can print the message. The catch aborts itself,
  // since other ranks may be waiting in a collective.
  alps::mpi::environment env(argc, argv);
#endif
  try {
    if (is_input_output_call(argc, argv)) {
      alps::solvers::cthyb_segment(alps::read_params_file(argv[1]), argv[2]);
      return 0;
    }
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
