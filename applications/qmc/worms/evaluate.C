/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2002-2003 by Matthias Troyer <troyer@comp-phys.org>,
*                            Simon Trebst <trebst@comp-phys.org>
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

void evaluate(const boost::filesystem::path& p, std::ostream& out, const bool write_xml) {
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
  alps::RealObsevaluator n  = sim.get_measurements()["Density"];
  alps::RealObsevaluator n2 = sim.get_measurements()["Density^2"];

  alps::RealObsevaluator kappa= beta*(n2 - numsites*n*n);  // add factor of beta
  kappa.rename("Compressibility");

  bool use_1D_stiffness = parms.defined("USE_1D_STIFFNESS") ? static_cast<bool>(parms["USE_1D_STIFFNESS"]) : false ; 
  //bool use_1D_stiffness = parms.defined("USE_1D_STIFFNESS") ? true : false ; 

  if (use_1D_stiffness){
    // winding number statistics
    alps::RealVectorObsevaluator wz = sim.get_measurements()["Winding number histogram"] ;
    alps::RealObsevaluator wz_m1(wz.slice(0)) ;
    alps::RealObsevaluator wz_0(wz.slice(1)) ;
    alps::RealObsevaluator wz_p1(wz.slice(2)) ;

    alps::RealObsevaluator rho_s_01=log(2.*wz_0/( wz_m1+wz_p1 )) ;
    rho_s_01 *= 2.*beta/L ;
    rho_s_01 = 1./rho_s_01 ; 
    rho_s_01.rename("Superfluid stiffness (1D estimator)") ;
    out<<rho_s_01<<"\n" ; 
    sim << rho_s_01 ; 
  }



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
    evaluate(p,std::cout,write_xml);
   }

/*  boost::filesystem::path p(argv[1]);
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
