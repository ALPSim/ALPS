/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2002-2008 by Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/scheduler.h>
#include <alps/alea/observables.hpp>
#include <fstream>
#include <filesystem>

void evaluate(const std::filesystem::path& p, const bool write_xml) {
  alps::ProcessList nowhere;
  alps::scheduler::MCSimulation sim(nowhere,p);

  // read in parameters
  alps::Parameters parms=sim.get_parameters();
  double beta=parms.defined("beta") ? alps::evaluate<double>("beta",parms) : 
         1./alps::evaluate<double>("T",parms);             
  alps::graph_helper<> graph(parms);
  double numsites = graph.num_sites();

  // determine specific heat
  const alps::alea::result_set measurements = sim.get_measurements();
  alps::alea::observable_result n  = measurements["n"];
  alps::alea::observable_result n2 = measurements["n^2"];
  alps::alea::observable_result cv = (n2-n*n-n)/numsites;
  std::cout << "Specific Heat: " << cv << "\n";
  sim.add_result("Specific Heat", cv);

  if (measurements.has("Density")) {
    // determine compressibility
    alps::alea::observable_result density  = measurements["Density"];
    alps::alea::observable_result density2 = measurements["Density^2"];

    alps::alea::observable_result kappa = beta*(density2 - numsites*density*density);
    std::cout << "Compressibility: " << kappa << "\n";
    sim.add_result("Compressibility", kappa);

    if(parms.value_or_default("MEASURE[Local Compressibility]",false)) {
      alps::alea::observable_result nl  = measurements["Local Density"];
      alps::alea::observable_result nlg  = measurements["Local Density * Global Density"];

      alps::alea::observable_result kappa_local = beta*numsites*(nlg - nl*density);
      std::cout << "Local Compressibility: " << kappa_local << "\n";
      sim.add_result("Local Compressibility", kappa_local);
    }

    if(parms.value_or_default("MEASURE[Site Compressibility]",false)) {
      alps::alea::observable_result nl  = measurements["Local Density"];
      alps::alea::observable_result nlc  = measurements["Integrated Local Density Correlations"];

      alps::alea::observable_result kappa_site = beta*(nl - nl*nlc);
      std::cout << "Site Compressibility: " << kappa_site << "\n";
      sim.add_result("Site Compressibility", kappa_site);
    }
  }
  // save
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
    std::filesystem::path p =
      std::filesystem::absolute(std::filesystem::path(argv[i]));
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
