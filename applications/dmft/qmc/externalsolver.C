 /*****************************************************************************
 *
 * ALPS DMFT Project
 *
 * Copyright (C) 2005 - 2009 by Emanuel Gull <gull@phys.columbia.edu>
 *                              Philipp Werner <werner@itp.phys.ethz.ch>,
 *                              Sebastian Fuchs <fuchs@theorie.physik.uni-goettingen.de>
 *                              Matthias Troyer <troyer@comp-phys.org>
 *
 *
 * THIS SOFTWARE NEEDS AN APPROPRIATE LICENSE BLOCK HERE
 *****************************************************************************/

/* $Id: externalsolver.C 360 2009-06-01 02:32:00Z gullc $ */

/// @file externalsolver.C
/// @brief implements the external solver
/// @sa ExternalSolver

#include "externalsolver.h"
#include "xml.h"
#include <alps/xml.h>
#include <alps/parser/parser.h>
#include <alps/utility/vectorio.hpp>
#include <cstdlib>
#include <fstream>
#include <boost/tuple/tuple.hpp>

#ifdef BOOST_MSVC
#include <io.h>
#endif

ImpuritySolver::result_type ExternalSolver::solve(
      const itime_green_function_t& G0
    , const alps::Parameters& parms) 
{
  char basename_in[]="alps_external_solver_XXXXXX";
  char basename_out[]="alps_external_solver_out_XXXXXX";
  mkstemp(basename_in); // the input file name
  mkstemp(basename_out); // the output file name
  std::string infile(basename_in), outfile(basename_out);
  if(parms.defined("TMPNAME")){
    infile=basename_in+std::string(".in.xml");
    outfile=basename_out+std::string(".out.xml");
  }
  boost::filesystem::path inpath(infile,boost::filesystem::native);
  boost::filesystem::path outpath(outfile,boost::filesystem::native);
 
  // write input XML file
  {
    alps::oxstream input(boost::filesystem::complete(inpath));
  
    input << alps::start_tag("SIMULATION") << parms;
	input << alps::start_tag("G0");
    write_itime(input,G0);
    input << alps::end_tag("G0");
    input << alps::end_tag("SIMULATION");
	
    // scope to force closing of file
  }
  call(infile,outfile);

  // read the output
  unsigned int N=(unsigned int)parms["N"];
  unsigned int sites     =(unsigned int)parms.value_or_default("SITES", 1);
  unsigned int flavors   =(unsigned int)parms.value_or_default("FLAVORS", 2);
 
  itime_green_function_t g(N+1, sites, flavors);
  
  {
    std::ifstream output(outfile.c_str());
    alps::XMLTag tag=alps::parse_tag(output); // skip outermost element <SIMULATION>
	
	// search for up and down spin Green's function
	tag=alps::parse_tag(output);
	while (tag.name != "GREENFUNCTION") {
	  if (tag.name=="/SIMULATION")
	    boost::throw_exception(std::runtime_error("Element <GREENFUNCTION> missing in output file"));
	  else
	    alps::skip_element(output,tag);
	  tag=alps::parse_tag(output);
	}

    read_itime(output,g);	
    
  } // scope to force closing of file
  
  boost::filesystem::remove(boost::filesystem::complete(outpath));
  
  return g;
}

MatsubaraImpuritySolver::result_type ExternalSolver::solve_omega(
              const matsubara_green_function_t& G0_omega
            , const alps::Parameters& parms)
{
  char basename_in[]="alps_external_solver_XXXXXX";
  char basename_out[]="alps_external_solver_out_XXXXXX";
  mkstemp(basename_in); // the input file name
  mkstemp(basename_out); // the output file name
  std::string infile(basename_in), outfile(basename_out);
  if(parms.defined("TMPNAME")){
    infile=basename_in+std::string(".in.xml");
    outfile=basename_out+std::string(".out.xml");
  }
  boost::filesystem::path inpath(infile,boost::filesystem::native);
  boost::filesystem::path outpath(outfile,boost::filesystem::native);
 
  // write input XML file
  {
    alps::oxstream input(boost::filesystem::complete(inpath));
  
    input << alps::start_tag("SIMULATION") << parms;
	input << alps::start_tag("G0");
    write_freq(input,G0_omega);
    input << alps::end_tag("G0");
    input << alps::end_tag("SIMULATION");
	
    // scope to force closing of file
  }
  
  call(infile,outfile);

  // read the output
  
  unsigned int n_matsubara=(unsigned int)parms["NMATSUBARA"];
  unsigned int n_tau=(unsigned int)parms["N"];
  unsigned int n_site     =(unsigned int)parms.value_or_default("SITES", 1);
  unsigned int n_orbital  =(unsigned int)parms.value_or_default("FLAVORS", 2);
  matsubara_green_function_t G_omega(n_matsubara, n_site, n_orbital);
  itime_green_function_t G_tau(n_tau+1, n_site, n_orbital);
  {
    std::ifstream output(outfile.c_str());
    alps::XMLTag tag=alps::parse_tag(output); // skip outermost element <SIMULATION>
    // search for up and down spin Green's function
    tag=alps::parse_tag(output);
    while (tag.name != "GREENFUNCTION_OMEGA") {
      if (tag.name=="/SIMULATION")
        boost::throw_exception(std::runtime_error("Element <GREENFUNCTION_OMEGA> missing in output file"));
      else
        alps::skip_element(output,tag);
      tag=alps::parse_tag(output);
    }
    read_freq(output,G_omega);	
  } // scope to force closing of file
  {
    std::ifstream output(outfile.c_str());
    alps::XMLTag tag=alps::parse_tag(output); // skip outermost element <SIMULATION>
	
    // search for up and down spin Green's function
    tag=alps::parse_tag(output);
    while (tag.name != "GREENFUNCTION_ITIME") {
      if (tag.name=="/SIMULATION")
        boost::throw_exception(std::runtime_error("Element <GREENFUNCTION_ITIME> missing in output file"));
      else
        alps::skip_element(output,tag);
      tag=alps::parse_tag(output);
    }
    read_itime(output,G_tau);	
  } // scope to force closing of file
  
  
  boost::filesystem::remove(boost::filesystem::complete(outpath));
  
  return std::make_pair(G_omega, G_tau);
}



void ExternalSolver::call(std::string const& infile, std::string const& outfile)
{
  boost::filesystem::path inpath(infile,boost::filesystem::native);
  boost::filesystem::path outpath(outfile,boost::filesystem::native);
  
  // call the external solver program
  std::string command = exe_.native_file_string() + " " + infile + " " + outfile;
  std::cerr << "Calling external solver " << exe_.native_file_string() << "\n\n";
  std::cerr <<" command is: "<<command<<std::endl;
  int result = std::system(command.c_str());
  if (result)
    boost::throw_exception(std::runtime_error("System error code " +boost::lexical_cast<std::string>(result) + " encountered when executing command:\n"+command));
  //std::cerr << "\nFinished call to external solver.\n";
	                                           
  boost::filesystem::remove(boost::filesystem::complete(inpath));

  if (!boost::filesystem::exists(outpath))
    boost::throw_exception(std::runtime_error("The external impurity solver failed to write the output file named " + outfile));

}
