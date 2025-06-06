/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2002-2008 by Matthias Troyer <troyer@comp-phys.org>
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

#include <alps/scheduler.h>
#include <alps/alea.h>
#include <fstream>
#include <boost/filesystem/operations.hpp>

void evaluate(const boost::filesystem::path& p, const bool write_xml) {
  alps::ProcessList nowhere;
  alps::scheduler::MCSimulation sim(nowhere,p);

  // read in parameters
  alps::Parameters parms=sim.get_parameters();
  double beta=parms.defined("beta") ? alps::evaluate<double>("beta",parms) : 
         1./alps::evaluate<double>("T",parms);             
  alps::graph_helper<> graph(parms);
  double numsites = graph.num_sites();

  // determine specific heat
  
  alps::RealObsevaluator n  = sim.get_measurements()["n"];
  alps::RealObsevaluator n2 = sim.get_measurements()["n^2"];
  alps::RealObsevaluator cv = (n2-n*n-n)/numsites;
  cv.rename("Specific Heat");
  std::cout << cv << "\n";
//  for (int i=0;i<cv.bin_number();++i)
//    std::cout << std::setprecision(20) << cv.bin_value(i)/cv.bin_size() << "\n";
	
  sim << cv;

  if (sim.get_measurements().has("Density")) {
	  // determine compressibility
	  alps::RealObsevaluator density  = sim.get_measurements()["Density"];
	  alps::RealObsevaluator density2 = sim.get_measurements()["Density^2"];

	  alps::RealObsevaluator kappa= (density2 - numsites*density*density);  // add factor of beta
	  kappa*=beta;
	  kappa.rename("Compressibility");
	  std::cout << kappa << "\n";
	  sim << kappa;

	  if(parms.value_or_default("MEASURE[Local Compressibility]",false)) {
		alps::RealVectorObsevaluator nl  = sim.get_measurements()["Local Density"];
		alps::RealVectorObsevaluator nlg  = sim.get_measurements()["Local Density * Global Density"];

		alps::RealVectorObsevaluator kappa_local= (nlg - nl*density);  // add factor of beta
		kappa_local *= beta*numsites;
		kappa_local.rename("Local Compressibility");
		std::cout << kappa_local << "\n";
		sim << kappa_local;
	  }

	  if(parms.value_or_default("MEASURE[Site Compressibility]",false)) {
		alps::RealVectorObsevaluator nl  = sim.get_measurements()["Local Density"];
		alps::RealVectorObsevaluator nlc  = sim.get_measurements()["Integrated Local Density Correlations"];

		alps::RealVectorObsevaluator kappa_site = nl - nl*nlc;  // add factor of beta
		kappa_site *= beta;
		kappa_site.rename("Site Compressibility");
		std::cout << kappa_site << "\n";
		sim << kappa_site;
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
