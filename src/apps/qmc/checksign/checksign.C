/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
    tag=alps::parse_tag(xml,true);
    while (tag.name!="/JOB") {
      if (tag.name=="TASK") {
        tag=alps::parse_tag(xml,true);
        while (tag.name !="/TASK") {
          if (tag.name=="INPUT") {
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
