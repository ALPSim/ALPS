/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2002-2003 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/cli.hpp>
#include <alps/utilities/copyright.hpp>
#include <alps/scheduler.h>
#include <alps/alea/observables.hpp>
#include <fstream>
#include <filesystem>

void evaluate(const std::filesystem::path& p, std::ostream& out, const bool write_xml) {
  alps::ProcessList nowhere;
  alps::scheduler::MCSimulation sim(nowhere,p);

  // read in parameters
  alps::Parameters parms=sim.get_parameters();
  double beta=parms.defined("beta") ? alps::evaluate<double>("beta",parms) : 
         1./alps::evaluate<double>("T",parms);     

  int L=alps::evaluate<double>("L",parms) ;
         
  alps::graph_helper<> graph(parms);
  double numsites = graph.num_sites();

  // determine compressibility
  const alps::alea::result_set measurements = sim.get_measurements();
  alps::alea::observable_result n  = measurements["Density"];
  alps::alea::observable_result n2 = measurements["Density^2"];

  alps::alea::observable_result kappa= beta*(n2 - numsites*n*n);  // add factor of beta

  bool use_1D_stiffness = parms.defined("USE_1D_STIFFNESS") ? static_cast<bool>(parms["USE_1D_STIFFNESS"]) : false ; 

  if (use_1D_stiffness){
    // winding number statistics: how often the winding number was -1, 0 and +1
    alps::alea::observable_result rho_s_01 = alps::alea::transform(
      [beta, L](const std::vector<double>& wz) {
        return std::vector<double>(1, 1./(2.*beta/L*std::log(2.*wz[1]/(wz[0]+wz[2]))));
      }, {measurements["Winding number histogram"]});
    out << "Superfluid stiffness (1D estimator): " << rho_s_01 << "\n";
    sim.add_result("Superfluid stiffness (1D estimator)", rho_s_01);
  }

  // output
  out << "Compressibility: " << kappa << "\n";
  sim.add_result("Compressibility", kappa);
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
  if (alps::handle_cli_information(argc, argv, "worm", [&] {
    std::cout << "Usage: " << argv[0] << " [--write-xml] inputfile [inputfile ...]\n";
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
  if (alps::cli_is_master()) alps::print_copyright(std::cout, "worm");

  for(; i<argc; i++)
   {
    std::filesystem::path p =
      std::filesystem::absolute(std::filesystem::path(argv[i]));
    evaluate(p,std::cout,write_xml);
   }

/*  std::filesystem::path p(argv[1]);
  std::string name=argv[1];
  if (argc==2)
    evaluate(p,std::cout);
  else {
    std::ofstream output(argv[2]);
    evaluate(p,output);
  }*/

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  std::exit(-5);
}
#endif
}
