/*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2010 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
 *
 *****************************************************************************/

/* $Id: solver_main.C 285 2008-01-03 15:34:39Z gullc $ */

#include "hirschfyesim.h"
#include <alps/parameter.h>
#include <alps/parameter_input.hpp>
#include <alps/utilities/copyright.hpp>
#include <alps/cli.hpp>
#include <alps/cli/options.hpp>
#include <alps/utilities/vectorio.hpp>
#include <boost/throw_exception.hpp>

bool parse_options(int argc, char** argv, std::string& infile, std::string& outfile)
{

  alps::cli::options desc("Allowed options");
  desc.flag("help", "produce help message")
    .flag("license,l", "print license conditions")
    .value("input-file", infile, "input file")
    .value("output-file", outfile, "output file")
    .positional("input-file").positional("output-file");
  if (alps::handle_cli_information(argc, argv, "hirschfye", [&] { std::cout << desc << "\n"; })) return false;
  desc.parse(argc, argv);

  if (infile.empty() || outfile.empty())
    throw std::invalid_argument("Expected input and output files");
  return true;
}

/// @brief The main program of the impurity solver
///
/// The program must be called with at least two command line parameters: the name of the input and output files.
/// Additional command line options are --help to print the usage information and --license to print license information

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
	try {
#endif
		std::string infile;
		std::string outfile;
		if (!parse_options(argc,argv,infile,outfile))
			return 0;
    alps::cli_mpi_guard mpi(argc, argv);
		// read parameters and G0
		
    alps::Parameters parms = alps::read_parameter_file(infile);
    parms["INFILE"]=infile;
    parms["OUTFILE"]=outfile;
    int N=(int)parms["NMATSUBARA"];
    int sites=parms.value_or_default("SITES", 1);
    int flavors=parms.value_or_default("FLAVORS", 2);
    
    alps::hdf5::archive ar(infile, "r");
    matsubara_green_function_t g0(N, sites, flavors); g0.read_hdf5(ar, "/G0");
    if (alps::cli_is_master()) {
      std::cout << "ALPS Hirsch-Fye solver for the single site impurity problem.\n\n";
      alps::print_copyright(std::cout, "hirschfye");
    }
    alps::scheduler::BasicFactory<HirschFyeSim,HirschFyeRun> factory;
    alps::ImpuritySolver solver(factory,parms,argc,argv);
		
		// write g into output file
    solver.solve_omega(g0,parms);
#ifndef BOOST_NO_EXCEPTIONS
	}
	catch (std::exception& exc) {
		std::cerr << exc.what() << "\n";
		return -1;
  }
	catch (...) {
		std::cerr << "Fatal Error: Unknown Exception!\n";
		return -2;
	}
#endif  
	return 0;
}
