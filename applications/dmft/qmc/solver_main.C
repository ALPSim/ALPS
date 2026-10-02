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
#include <alps/utility/copyright.hpp>
#include <alps/utility/cli.hpp>
#include <alps/utility/vectorio.hpp>
#include <boost/throw_exception.hpp>
#include <boost/program_options.hpp>
#ifdef BOOST_MSVC
#include <direct.h>
#endif

bool parse_options(int argc, char** argv, std::string& infile, std::string& outfile)
{

	namespace po = boost::program_options;
	
	po::options_description desc("Allowed options");
	desc.add_options()
  ("help", "produce help message")
  ("license,l", "print license conditions") 
  ("input-file", po::value<std::string>(&infile), "input file")
  ("output-file", po::value<std::string>(&outfile), "output file");
  if (alps::handle_cli_information(argc, argv, "hirschfye", [&] { std::cout << desc << "\n"; })) return false;
	po::positional_options_description p;
	p.add("input-file", 1);
	p.add("output-file", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);    
	
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
		
    alps::hdf5::archive ar(infile, "r");
		alps::Parameters parms;
    ar["/parameters"] >> parms;
    parms["INFILE"]=infile;
    parms["OUTFILE"]=outfile;
    int N=(int)parms["NMATSUBARA"];
    int sites=parms.value_or_default("SITES", 1);
    int flavors=parms.value_or_default("FLAVORS", 2);
    
    matsubara_green_function_t g0(N, sites, flavors); g0.read_hdf5(ar, "/G0");
    if (alps::cli_is_master()) {
      std::cout << "ALPS Hirsch-Fye solver for the single site impurity problem.\n\n";
      alps::print_copyright(std::cout, "hirschfye");
    }
    alps::scheduler::BasicFactory<HirschFyeSim,HirschFyeRun> factory;
    alps::ImpuritySolver solver(factory,argc,argv,true);
		
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
