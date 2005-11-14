/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005 by Philipp Werner <werner@comp-phys.org>,
*                       Matthias Troyer <troyer@comp-phys.org>
*                       Emanuel Gull <gullc@comp-phys.org>
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
#include <alps/xml.h>
#include <alps/parser/parser.h>
#include <alps/vectorio.h>
#include <cstdlib>
#include <fstream>

multiple_vector_type ExternalSolver::solve(const multiple_vector_type& G0, const alps::Parameters& parms) 
{
  std::string infile = std::tmpnam(0); // the input file name
  std::string outfile = std::tmpnam(0); // the output file name
 
  boost::filesystem::path inpath(infile,boost::filesystem::native);
  boost::filesystem::path outpath(outfile,boost::filesystem::native);
 
  // write input XML file
  {
    alps::oxstream input(inpath);
  
    input << alps::start_tag("SIMULATION") << parms;
	input << alps::start_tag("G0") << alps::attribute("spin","up");
	for (std::size_t i=0;i<G0.first.size();++i)
	  input << G0.first[i];
    input << alps::end_tag("G0");
	input << alps::start_tag("G0") << alps::attribute("spin","down");
	for (std::size_t i=0;i<G0.second.size();++i)
	  input << G0.second[i];
    input << alps::end_tag("G0");
		  
    input << alps::end_tag("SIMULATION");
	
    // scope to force closing of file
  }
  
  // call the external solver program
  std::string command = exe_.native_file_string() + " " + infile + " " + outfile;
  std::cerr << "Calling external solver " << exe_.native_file_string() << "\n\n";
  int result = std::system(command.c_str());
  if (result)
    boost::throw_exception(std::runtime_error("System error code " +boost::lexical_cast<std::string>(result) + " encountered when executing command:\n"+command));
  std::cerr << "\nFinished call to external solver.\n";
	                                           
  boost::filesystem::remove(inpath);

  // read the Green's function g from the output file
  multiple_vector_type g;
  
  if (!boost::filesystem::exists(outpath))
    boost::throw_exception(std::runtime_error("The external impurity solver failed to write the output file named " + outfile));

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

	bool first_is_up = (tag.attributes["spin"]=="up");
	if (tag.type==alps::XMLTag::OPENING) {
	  alps::read_vector_resize(alps::parse_content(output), first_is_up ? g.first : g.second);
      alps::skip_element(output,tag);
	}
	
    // search for second
	tag=alps::parse_tag(output);
	while (tag.name != "GREENFUNCTION") {
	  if (tag.name=="/SIMULATION")
	    boost::throw_exception(std::runtime_error("Element <GREENFUNCTION> missing in output file"));
	  else
	    alps::skip_element(output,tag);
	  tag=alps::parse_tag(output);
	}
	if (tag.type==alps::XMLTag::OPENING) {
	  alps::read_vector_resize(alps::parse_content(output), first_is_up ? g.second : g.first);
      alps::skip_element(output,tag);
	}
	
  } // scope to force closing of file
  
  boost::filesystem::remove(outpath);
  
  return g;
}
