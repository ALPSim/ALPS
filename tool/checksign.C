/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

#include <alps/model.h>
#include <alps/lattice.h>
#include <alps/parameterlist.h>
#include <alps/xml.h>
#include <iostream>
#include <string>

void check_parameters(const std::string& n, int i, alps::Parameters& p)
{
  alps::ModelLibrary models(p);
  alps::graph_factory<> lattices(p);
  alps::HamiltonianDescriptor<short> ham(models.hamiltonian(p["MODEL"]));
  p.copy_undefined(ham.default_parameters());
  ham.set_parameters(p);
  std::cout << n;
  if (i)
    std::cout << ", task " << i;
    
  if (has_sign_problem(ham,lattices.graph(),models.simple_operators(),p))
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

void check_xmlfile(int i, const std::string& name)
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
            check_xmlfile(j,tag.attributes["file"]);
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

  if (argc<2) {
    std::cerr << "Usage: " << argv[0] << " inputfile [inputfile ...]]\n";
    std::exit(-1);
  }
  for (int i=1;i<argc;++i) {    
    char c;
    { // check for XML file
      std::ifstream in(argv[i]);
      in >> c;
    }

    if (c=='<') // we have an XML file
      check_xmlfile(1,argv[i]);
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
