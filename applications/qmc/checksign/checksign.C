/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>
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

/* $Id$ */

#include <alps/model.h>
#include <alps/lattice.h>
#include <alps/utility/copyright.hpp>
#include <alps/xml.h>
#include <cstring>
#include <iostream>
#include <string>

void check_parameters(const std::string& n, int i, alps::Parameters& p)
{
  alps::graph_helper<> lattice(p);
  alps::model_helper<> models(lattice, p);
  std::cout << n;
  if (i)
    std::cout << ", task " << i;
    
  if (alps::has_sign_problem(models.model(),lattice,p))
    std::cout << ": SIGN PROBLEM\n";
  else
    std::cout << ": OK\n";
}

void check_parameterfile(const std::string& name)
{
  alps::ParameterList parms;
  { // scope for ifstream lifetime
    std::ifstream in(name.c_str());
    in >> parms;
  }
  if (parms.size()==0) { // we got a single simulation parameter set
    std::ifstream in(name.c_str());
    alps::Parameters p;
    in >> p;
    check_parameters(name,0,p);
  }
  else
    for (int i=0;i<parms.size();++i)
      check_parameters(name,i+1,parms[i]);
}

void check_xmlfile(const std::string& name)
{
  alps::Parameters p;
  std::ifstream xml(name.c_str());
  alps::XMLTag tag=alps::parse_tag(xml,true);
  if (tag.name=="JOB") {
    int j=0;
    tag=alps::parse_tag(xml,true);
    while (tag.name!="/JOB") {
      if (tag.name=="TASK") {
        tag=alps::parse_tag(xml,true);
        while (tag.name !="/TASK") {
          if (tag.name=="INPUT") {
            ++j;         
            check_xmlfile(tag.attributes["file"]);
          }
          else
            alps::skip_element(xml,tag);
          tag=alps::parse_tag(xml,true);
        }
      }
      else
        alps::skip_element(xml,tag);
      tag=alps::parse_tag(xml,true);
    }
  }
  else if (tag.name=="SIMULATION") {
    tag=alps::parse_tag(xml,true);
    while (tag.name!="/SIMULATION") {
      if (tag.name=="PARAMETERS") {
        // read parameters
        alps::Parameters p;
        p.read_xml(tag,xml);
        // check
        check_parameters(name,0,p);
      }
      else
        alps::skip_element(xml,tag);
      tag=alps::parse_tag(xml,true);
    }
  }
  else
    boost::throw_exception(std::runtime_error("Got unknown XML file starting with element " + tag.name));  
}


int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

  std::cout << "ALPS application to check for a sign problem in a quantum model\n"
            << "  available from http://alps.comp-phys.org/\n"
            << "  copyright (c) 2003-2007 by Matthias Troyer <troyer@comp-phys.org>\n\n";
  alps::print_copyright(std::cout);
  
  if (argc<2) {
    std::cerr << "Usage: " << argv[0] << " [-l] inputfile [inputfile ...]]\n";
    std::exit(-1);
  }
  for (int i=1;i<argc;++i) {   
    if (!std::strcmp(argv[i],"-l")) {
      alps::print_license(std::cout);
      continue;
    }
    char c;
    { // check for XML file
      std::ifstream in(argv[i]);
      in >> c;
    }

    if (c=='<') // we have an XML file
      check_xmlfile(argv[i]);
    else // we have a text file
      check_parameterfile(argv[i]);
  }
    

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
}
