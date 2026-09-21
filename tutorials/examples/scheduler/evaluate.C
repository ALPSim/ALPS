/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/scheduler/montecarlo.h>
#include <alps/alea.h>
#include <fstream>

void evaluate(const boost::filesystem::path& p, std::ostream& out) {
  alps::ProcessList nowhere;
  alps::scheduler::MCSimulation sim(nowhere,p);
  alps::RealObsevaluator m2=sim.get_measurements()["Magnetization^2"];
  alps::RealObsevaluator m4=sim.get_measurements()["Magnetization^4"];
#ifdef ALPS_HAVE_VALARRAY
  alps::RealVectorObsevaluator corr=sim.get_measurements()["Correlations"];
#endif
  alps::RealObsevaluator binder=m4/(m2*m2);
 
  binder.rename("Binder cumulant of Magnetization");
#ifdef ALPS_HAVE_VALARRAY
  out << corr << "\n";
#endif
  out << binder << "\n";
  sim.addObservable(binder);
  sim.checkpoint(p);
}

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::scheduler::SimpleMCFactory<alps::scheduler::DummyMCRun> factory;
  alps::scheduler::init(factory);
  
  if (argc<2 || argc>3) {
    std::cerr << "Usage: " << argv[0] << " inputfile [outputbasename]\n";
    std::exit(-1);
  }
  boost::filesystem::path p(argv[1]);
  std::string name=argv[1];
  if (argc==2)
    evaluate(p,std::cout);
  else {
    std::ofstream output(argv[2]);
    evaluate(p,output);
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
