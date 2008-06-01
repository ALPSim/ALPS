/*******************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005-2006 by Philipp Werner <werner@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*                            Emanuel Gull <gullc@comp-phys.org>
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

/* $Id$ */

/// @file externalsolver.C
/// @brief implements the external solver
/// @sa ExternalSolver

#include "externalsolver.h"
#include "xml.h"
#include <alps/xml.h>
#include <alps/parser/parser.h>
#include <alps/vectorio.h>
#include <cstdlib>
#include <fstream>
#include <boost/tuple/tuple.hpp>

ImpuritySolver::result_type ExternalSolver::solve(
      const itime_green_function_t& G0
    , const alps::Parameters& parms) 
{
  std::string infile = std::tmpnam(0); // the input file name
  std::string outfile = std::tmpnam(0); // the output file name
  if(parms.defined("TMPNAME")){
    infile=parms["TMPNAME"]+std::string(".in.xml");
    outfile=parms["TMPNAME"]+std::string(".out.xml");
  }
  boost::filesystem::path inpath(infile,boost::filesystem::native);
  boost::filesystem::path outpath(outfile,boost::filesystem::native);
 
  // write input XML file
  {
    alps::oxstream input(inpath);
  
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
  
  boost::filesystem::remove(outpath);
  
  return g;
}

MatsubaraImpuritySolver::result_type ExternalSolver::solve_omega(
              const matsubara_green_function_t& G0_omega
            , const alps::Parameters& parms)
{
  std::string infile = std::tmpnam(0); // the input file name
  std::string outfile = std::tmpnam(0); // the output file name
  if(parms.defined("TMPNAME")){
    infile=parms["TMPNAME"]+std::string(".in.xml");
    outfile=parms["TMPNAME"]+std::string(".out.xml");
  }

  boost::filesystem::path inpath(infile,boost::filesystem::native);
  boost::filesystem::path outpath(outfile,boost::filesystem::native);
 
  // write input XML file
  {
    alps::oxstream input(inpath);
  
    input << alps::start_tag("SIMULATION") << parms;
	input << alps::start_tag("G0");
    write_freq(input,G0_omega);
    input << alps::end_tag("G0");
    input << alps::end_tag("SIMULATION");
	
    // scope to force closing of file
  }
  
  std::cout<<"finished writing input xml file"<<std::endl;
  std::cout<<"input file name was: "<<infile<<std::endl;
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
    std::cout<<"reading in freq"<<std::endl;
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
  
  
  boost::filesystem::remove(outpath);
  
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
	                                           
  boost::filesystem::remove(inpath);

  if (!boost::filesystem::exists(outpath))
    boost::throw_exception(std::runtime_error("The external impurity solver failed to write the output file named " + outfile));

}
