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
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

/* $Id: solver_main.C 285 2008-01-03 15:34:39Z gullc $ */

#include "hirschfyesim.h"
#include <alps/parameter.h>
#include <alps/utility/copyright.hpp>
#include <alps/utility/vectorio.hpp>
#include <boost/throw_exception.hpp>
#include <boost/program_options.hpp>
#ifdef BOOST_MSVC
#include <direct.h>
#endif

bool parse_options(int argc, char** argv, std::string& infile, std::string& outfile)
{
	std::cout << "ALPS Hirsch-Fye solver for the single site impurity problem.\n\n";
	alps::print_copyright(std::cout);
	
	namespace po = boost::program_options;
	
	po::options_description desc("Allowed options");
	desc.add_options()
  ("help", "produce help message")
  ("license,l", "print license conditions") 
  ("input-file", po::value<std::string>(&infile), "input file")
  ("output-file", po::value<std::string>(&outfile), "output file");
	po::positional_options_description p;
	p.add("input-file", 1);
	p.add("output-file", 1);
	
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
	po::notify(vm);    
	
	bool valid=true;
	
	if (vm.count("help")) {
		std::cout << desc << "\n";
		valid=false;
	}
	if (vm.count("license")) {
		alps::print_license(std::cout);
		valid=false;
	}
	return valid;
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
