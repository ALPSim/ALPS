/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2002-2003 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
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

#include <alps/scheduler.h>
#include <alps/alea.h>
#include <fstream>

void evaluate(const boost::filesystem::path& p, std::ostream& out) {
  alps::ProcessList nowhere;
  alps::scheduler::MCSimulation sim(nowhere,p);

  // read in parameters
  alps::Parameters parms=sim.get_parameters();
  double beta=parms.defined("beta") ? alps::evaluate<double>("beta",parms) : 
         1./alps::evaluate<double>("T",parms);     
         
  alps::graph_helper<> graph(parms);
  double numsites = graph.num_sites();

  // determine compressibility
  alps::RealObsevaluator n  = sim.get_measurements()["Particle number"];
  alps::RealObsevaluator n2 = sim.get_measurements()["Particle number^2"];

  alps::RealObsevaluator kappa= beta*(n2 - numsites*n*n);  // add factor of beta
  kappa.rename("Compressibility");

/*
  alps::RealVectorObsevaluator IntervalStatistics 
    = sim.get_measurements()["Statistics time intervals"];

  alps::RealObsevaluator intervals;
  alps::RealObsevaluator ratio = 0;
  for(int i=1; i<IntervalStatistics.value().size(); ++i) {
    intervals += IntervalStatistics.slice(i) * i;
    ratio += IntervalStatistics.slice(i) / i;
  }
*/

  // output
  out << kappa << "\n";
   //   << intervals << "\n"
   //   << ratio << "\n";
  sim << kappa;
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
  boost::filesystem::path p(argv[1],boost::filesystem::native);
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
