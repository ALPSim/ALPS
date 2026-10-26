/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2010 by Sebastian  Fuchs <fuchs@comp-phys.org>
*                       Thomas Pruschke <pruschke@comp-phys.org>
*                       Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/solvers.hpp>
#include <alps/cli.hpp>
#include <alps/utilities/copyright.hpp>
#include <alps/ngs/mcoptions.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/params.hpp>
#include <boost/lexical_cast.hpp>
#include <exception>
#include <iostream>

int main(int argc, char** argv) {
  try {
    alps::mcoptions options(argc, argv);
    if (!options.valid) return 0;
    alps::params parms;
    {
      alps::hdf5::archive input(options.input_file);
      input["/parameters"] >> parms;
    }
    parms.define<std::string>("BASENAME", options.output_file, "base name of the output file");
    std::string output_file = parms["BASENAME"].as<std::string>() + ".out.h5";
    alps::cli_mpi_guard mpi(argc, argv);
    if (alps::cli_is_master()) alps::print_copyright(std::cout);
    alps::solvers::maxent(parms, output_file);
    return 0;
  } catch (std::exception const& error) {
    std::cerr << "maxent: " << error.what() << '\n';
    return 1;
  }
}
