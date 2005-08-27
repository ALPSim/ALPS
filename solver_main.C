/*****************************************************************************
*
* ALPS DMFT Project
*
* Copyright (C) 2005 by Philipp Werner <werner@comp-phys.org>,
*                       Matthias Troyer <troyer@comp-phys.org>
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

#include "hirschfye.h"
#include <alps/parameters.h>
#include <alps/copyright.h>
#include <alps/vectorio.h>
#include <boost/throw_exception.hpp>
#include <boost/program_options.hpp>


void read_input(const std::string& inputfile, multiple_vector_type& G0,  alps::Parameters& parms)
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

	bool first_is_up = (tag.attributes["spin"]=="up");
	if (tag.type==alps::XMLTag::OPENING) {
	  alps::read_vector_resize(alps::parse_content(input), first_is_up ? G0.first : G0.second);
      alps::skip_element(input,tag);
	}
	
    // search for second
	tag=alps::parse_tag(input);
	while (tag.name != "G0") {
	  if (tag.name=="/SIMULATION")
	    boost::throw_exception(std::runtime_error("Element <G0> missing in input file"));
	  else
	    alps::skip_element(input,tag);
	  tag=alps::parse_tag(input);
	}
	
	if (tag.type==alps::XMLTag::OPENING) {
	  alps::read_vector_resize(alps::parse_content(input), first_is_up ? G0.second : G0.first);
      alps::skip_element(input,tag);
	}
}


void write_output(const std::string& outputfile, const multiple_vector_type& G0, const multiple_vector_type& G, const alps::Parameters& parms)
{
  alps::oxstream output(outputfile);
  
  output << alps::start_tag("SIMULATION") << parms;
  output << alps::start_tag("G0") << alps::attribute("spin","up");
  for (std::size_t i=0;i<G0.first.size();++i)
	output << G0.first[i];
  output << alps::end_tag("G0");
  output << alps::start_tag("G0") << alps::attribute("spin","down");
  for (std::size_t i=0;i<G0.second.size();++i)
	output << G0.second[i];
  output << alps::end_tag("G0");
  output << alps::start_tag("GREENFUNCTION") << alps::attribute("spin","up");
  for (std::size_t i=0;i<G.first.size();++i)
	output << G.first[i];
  output << alps::end_tag("GREENFUNCTION");
  output << alps::start_tag("GREENFUNCTION") << alps::attribute("spin","down");
  for (std::size_t i=0;i<G.second.size();++i)
	output << G.second[i];
  output << alps::end_tag("GREENFUNCTION");
  output << alps::end_tag("SIMULATION");
}


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
    ("output-file", po::value<std::string>(&outfile), "input file");
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
  
  multiple_vector_type G0;
  alps::Parameters parms;
  read_input(infile,G0,parms);
  
  HirschFyeSolver solver;

  // write g into output file
  multiple_vector_type G = solver.solve(G0,parms);
  write_output(outfile,G0,G,parms);
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
