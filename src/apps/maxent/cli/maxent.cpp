/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1998-2026 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/solvers.hpp>
#include <alps/cli.hpp>
#include <alps/ngs/mcoptions.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/params.hpp>
#include <boost/lexical_cast.hpp>
#include <exception>
#include <iostream>

int main(int argc, char** argv) {
  try {
    alps::mcoptions options(argc, argv, "maxent");
    if (!options.valid) return 0;
    alps::params parms;
    {
      alps::hdf5::archive input(options.input_file);
      input["/parameters"] >> parms;
    }
    parms.define<std::string>("BASENAME", options.output_file, "base name of the output file");
    std::string output_file = parms["BASENAME"].as<std::string>() + ".out.h5";
    alps::cli_mpi_guard mpi(argc, argv);
    alps::solvers::maxent(parms, output_file);
    return 0;
  } catch (std::exception const& error) {
    std::cerr << "maxent: " << error.what() << '\n';
    return 1;
  }
}
