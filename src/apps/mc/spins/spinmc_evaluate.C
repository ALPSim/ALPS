/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2003-2004 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/cli.hpp>
#include <alps/utilities/copyright.hpp>
#include <alps/config.h>
#include <alps/alea/observables.hpp>
#include <alps/scheduler.h>

#include <filesystem>

void evaluate(const std::filesystem::path& p, const bool write_xml) {
  typedef alps::alea::observable_result result;
  alps::ProcessList nowhere;
  alps::scheduler::MCSimulation sim(nowhere,p);
  const alps::alea::result_set m_in=sim.get_measurements();
  alps::Parameters parms=sim.get_parameters();
  alps::graph_helper<> graph(parms);
  double numsites=graph.num_sites();
  double beta=parms.defined("beta") ? alps::evaluate<double>("beta",parms) : 1./alps::evaluate<double>("T",parms);
  
  result obse_e = m_in["beta * Energy / sqrt(N)"];
  result obse_e2 = m_in["(beta * Energy)^2 / N"];
  sim.add_result("Specific Heat", obse_e2 - obse_e * obse_e);
  
  result m = m_in["|Magnetization|"];
  result m2 = m_in["Magnetization^2"];
  sim.add_result("Binder Cumulant U2", m2/(m*m));
  sim.add_result("Connected Susceptibility", beta*numsites*(m2-m*m));
  
  if (m_in.has("Magnetization^4")) {
    result m4 = m_in["Magnetization^4"];
    sim.add_result("Binder Cumulant", m4/(m2*m2));
  
    if (m_in.has("E.Magnetization^4")) {
      result e = m_in["Energy"];
      result em2 = m_in["E.Magnetization^2"];
      result em4 = m_in["E.Magnetization^4"];
      sim.add_result("Binder Cumulant slope", beta*beta*((em4-e*m4)/(m2*m2)-2*m4*(em2-e*m2)/(m2*m2*m2)));
      sim.add_result("Magnetization^2 slope", beta*beta*((em2-e*m2)));
      sim.add_result("Magnetization^4 slope", beta*beta*((em4-e*m4)));
    }
  }
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
  if (alps::handle_cli_information(argc, argv, "spinmc", [&] {
    std::cout << "Usage: " << argv[0] << " [--write-xml] inputfile [inputfile ...]\n"
              << "\"Binder Cumulant\" is the moment ratio <m^4>/<m^2>^2, \"Binder Cumulant U2\" is\n"
              << "<m^2>/<|m|>^2. \"Connected Susceptibility\" is beta N (<m^2> - <|m|>^2). spinmc's\n"
              << "\"Susceptibility\" is beta N <m^2>/n for n spin components.\n";
  })) return 0;
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


  if (i >= argc) throw std::invalid_argument("Expected an input file");
  alps::cli_mpi_guard mpi(argc, argv);
  if (alps::cli_is_master()) alps::print_copyright(std::cout, "spinmc");

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
