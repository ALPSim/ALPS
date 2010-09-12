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
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/* $Id: solver_main.C 285 2008-01-03 15:34:39Z gullc $ */

#include "hirschfyesim.h"
#include "xml.h"
#include <alps/parameter.h>
#include <alps/utility/copyright.hpp>
#include <alps/utility/vectorio.hpp>
#include <boost/throw_exception.hpp>
#include <boost/program_options.hpp>
#ifdef BOOST_MSVC
#include <direct.h>
#endif

/// @file solver_main.C
/// @brief example main program for a Hirsch-Fye solver


/// @brief parse the XML input file for an impurity solver
///
/// @param inputfile name of the input XML file
/// @param G0 contains on return the Weiss field read from the input file, in Matsubara frequencies
/// @param parms contains on return the simulation parameter read from the input file

matsubara_green_function_t read_input(const std::string& inputfile,  alps::Parameters& parms)
{
	std::ifstream input(inputfile.c_str());
	alps::XMLTag tag=alps::parse_tag(input); // skip outermost element <SIMULATION>
	
	// search for up and down spin Green's function
	tag=alps::parse_tag(input);
	while (tag.name != "G0") {
		if (tag.name=="PARAMETERS")
			parms.read_xml(tag,input);
		else if (tag.name=="/SIMULATION")
			boost::throw_exception(std::runtime_error("Element <G0> missing in input file"));
		else
			alps::skip_element(input,tag);
		tag=alps::parse_tag(input);
	}
        int N=(int)parms["NMATSUBARA"];
        int flavors=parms.value_or_default("FLAVORS", 2);
        int sites=parms.value_or_default("SITES", 1);
        matsubara_green_function_t G0(N, sites, flavors);
        read_freq(input,G0);
        return G0;
}


/// @brief write the XML output file for an impurity solver
///
/// @param outputfile name of the input XML file
/// @param G0 the Weiss field used in the solver
/// @param G the Green's function calculated by the impurity solver
/// @param parms the simulation parameter used in the solver
void write_output(const std::string& outputfile, const matsubara_green_function_t& G0, const matsubara_green_function_t& G_omega, const itime_green_function_t &G_tau, const alps::Parameters& parms)
{
	alps::oxstream output(outputfile);
	
	output << alps::start_tag("SIMULATION") << parms;
    
    output << alps::start_tag("G0");
    write_freq(output,G0);
    output << alps::end_tag("G0");
    std::cout<<"writing output. sizes are: "<<G_omega.nfreq()<<" "<<G_tau.ntime()<<std::endl;
    output << alps::start_tag("GREENFUNCTION_OMEGA") ;
    write_freq(output,G_omega);
    output << alps::end_tag("GREENFUNCTION_OMEGA");
    output << alps::start_tag("GREENFUNCTION_ITIME") ;
    write_itime(output,G_tau);
    output << alps::end_tag("GREENFUNCTION_ITIME");
	output << alps::end_tag("SIMULATION");
}

/// @brief parse the command line options
///
/// @param argc the number of arguments
/// @param argv the array of strings containing the arguments
/// @param infile contains on return the name of the input file, which was the first argument
/// @param outfile contains on return the name of the output file, which was the second argument
/// @return returns true if the simulation should continue to run and false otherwise. 
///         False is returned if the user just wants to see the help or license information
///         by specifying the --help or --license option

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
        // set working directory
        boost::filesystem::path p(outfile,boost::filesystem::native);
        chdir(p.branch_path().native_file_string().c_str());
		alps::Parameters parms;
		matsubara_green_function_t G0=read_input(infile,parms);
                alps::scheduler::BasicFactory<HirschFyeSim,HirschFyeRun> factory;
                alps::ImpuritySolver solver(factory,argc,argv);
		
		// write g into output file
                std::pair<matsubara_green_function_t, itime_green_function_t> G = solver.solve_omega(G0,parms);
		write_output(outfile,G0,G.first,G.second, parms);
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
