/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2004 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/config.h>
#include <alps/alea.h>
#include <alps/scheduler.h>

#include <boost/filesystem/operations.hpp>

void evaluate(const boost::filesystem::path& p, const bool write_xml) {
  alps::ProcessList nowhere;
  alps::scheduler::MCSimulation sim(nowhere,p);
  const alps::ObservableSet& m_in=sim.get_measurements();
  alps::Parameters parms=sim.get_parameters();
  alps::graph_helper<> graph(parms);
  double numsites=graph.num_sites();
  double beta=parms.defined("beta") ? alps::evaluate<double>("beta",parms) : 1./alps::evaluate<double>("T",parms);
  
  alps::RealObsevaluator obse_e = m_in["beta * Energy / sqrt(N)"];
  alps::RealObsevaluator obse_e2 = m_in["(beta * Energy)^2 / N"];
  alps::RealObsevaluator eval("Specific Heat");
  eval = (obse_e2 - obse_e * obse_e);
  sim << eval;
  
  alps::RealObsevaluator m = m_in["|Magnetization|"];
  alps::RealObsevaluator m2 = m_in["Magnetization^2"];
  alps::RealObsevaluator u2("Binder Cumulant U2");
  u2 = m2/(m*m);
  sim << u2;
  alps::RealObsevaluator chi("Connected Susceptibility");
  chi=beta*numsites*(m2-m*m);
  sim << chi;
  
  if (m_in.has("Magnetization^4")) {
  alps::RealObsevaluator m4 = m_in["Magnetization^4"];
  alps::RealObsevaluator u4("Binder Cumulant");
  u4 = m4/(m2*m2);
  sim << u4;
  
  if (m_in.has("E.Magnetization^4")) {
  alps::RealObsevaluator e = m_in["Energy"];
  alps::RealObsevaluator em2 = m_in["E.Magnetization^2"];
  alps::RealObsevaluator em4 = m_in["E.Magnetization^4"];
  alps::RealObsevaluator du4("Binder Cumulant slope");
  alps::RealObsevaluator dm4("Magnetization^4 slope");
  alps::RealObsevaluator dm2("Magnetization^2 slope");
	du4=beta*beta*((em4-e*m4)/(m2*m2)-2*m4*(em2-e*m2)/(m2*m2*m2));
	dm2=beta*beta*((em2-e*m2)); dm4=beta*beta*((em4-e*m4)); 
  sim << du4; sim << dm2; sim << dm4;
 } }
   sim.checkpoint(p,write_xml);
}

int main(int argc, char** argv)
{
  int i;
  char write_xml_flag[]="--write-xml";
  bool write_xml;
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif
  alps::scheduler::SimpleMCFactory<alps::scheduler::DummyMCRun> factory;
  alps::scheduler::init(factory);
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " [--write-xml] inputfile1 [intputfile2 [...]]\n";
    std::exit(-1);
  }

  if (strcmp(write_xml_flag,argv[1])==0)  {
   write_xml=true;
   i=2;
  }
  else {
   write_xml=false;
   i=1;
  }


  for(; i<argc; i++)
   {
    boost::filesystem::path p =
      boost::filesystem::absolute(boost::filesystem::path(argv[i]));
    evaluate(p,write_xml);
   }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  std::exit(-5);
}
#endif
}
