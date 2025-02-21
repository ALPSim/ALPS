/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2004-2006 by Stefan Wessel <wessel@comp-phys.org>
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

#include <alps/alea.h>
#include <alps/plot.h>
#include <alps/scheduler.h>
#include <boost/filesystem/operations.hpp>
#include <map>
#include <string>
#include <utility>
#include <valarray>

std::pair<const std::basic_string<char,std::char_traits<char>,std::allocator<char> >,unsigned long> _buggy;

void evaluate(const boost::filesystem::path& p, const alps::Parameters new_parms) {
  alps::ProcessList nowhere;
  alps::scheduler::MCSimulation sim(nowhere,p);
  alps::Parameters parms=sim.get_parameters();  
  alps::graph_helper<> lattice(parms);
  int Ns=num_sites(lattice.graph());
  double T_min=new_parms.value_or_default("T_MIN", parms.value_or_default("T_MIN",0.1));
  double T_max=new_parms.value_or_default("T_MAX", parms.value_or_default("T_MAX",10.));
  double T_delta=new_parms.value_or_default("DELTA_T", parms.value_or_default("DELTA_T",0.1));
  alps::RealVectorObsevaluator g_=sim.get_measurements()["Coefficients"];
  std::valarray<double> g=g_.mean();
  alps::RealObsevaluator offset=sim.get_measurements()["Offset"];
  alps::plot::Set<double> energy_set;
  alps::plot::Set<double> free_energy_set;
  alps::plot::Set<double> entropy_set;
  alps::plot::Set<double> specific_heat_set;

// calculate thermdynamic properties

  for (double T=T_min;T<T_max+T_delta/2.;T+=T_delta) {    
    double logbeta=-log(T);
    double maxx=g[0];
    for (int i=1;i<g.size();++i)  // avoid overflow
      maxx=std::max(g[i]+logbeta*i,maxx);
    double Z=0;
    double Sn=0;
    double Sn2=0;
    for (int i=0;i<g.size();++i) {
      double factor=exp(g[i]+logbeta*i-maxx);
      Z+=factor;
      Sn+=i*factor;
      Sn2+=i*i*factor;
    };
    double energy       =(offset.mean()-T*Sn/Z)/Ns;
    double free_energy  =offset.mean()/Ns-T*(log(Z)+maxx)/Ns;
    double entropy      =(energy-free_energy)/T;
    double specific_heat=(Sn2/Z-(Sn/Z)*(Sn/Z)-Sn/Z)/Ns;
    energy_set        << T << energy;
    free_energy_set   << T << free_energy;
    entropy_set       << T << entropy;
    specific_heat_set << T << specific_heat;
  };
   std::string ss=p.string();

  ss.erase(ss.rfind(".out.xml"),8);
  if (1) {
    alps::plot::Plot<double> my_plot("Energy Density versus Temperature",parms);
    my_plot.set_labels("Temperature","Energy Density");
    my_plot << energy_set;
    alps::oxstream my_ox(ss+".plot.energy.xml");
    my_ox << my_plot;
  }
  if (1) {
    alps::plot::Plot<double> my_plot("Free Energy Density versus Temperature",parms);
    my_plot.set_labels("Temperature","Free Energy Density");
    my_plot << free_energy_set;
    alps::oxstream my_ox(ss+".plot.free_energy.xml");
    my_ox << my_plot;
  }
  if (1) {
    alps::plot::Plot<double> my_plot("Entropy Density versus Temperature",parms);
    my_plot.set_labels("Temperature","Entropy Density");
    my_plot << entropy_set;
    alps::oxstream my_ox(ss+".plot.entropy.xml");
    my_ox << my_plot;
  }
  if (1) {
    alps::plot::Plot<double> my_plot("Specific Heat per Site versus Temperature",parms);
    my_plot.set_labels("Temperature","Specific Heat per Site");
    my_plot << specific_heat_set;
    alps::oxstream my_ox(ss+".plot.specific_heat.xml");
    my_ox << my_plot;
  }
 
  // now calcuate magnetic  observables

  if (bool(parms.value_or_default("MEASURE_MAGNETIC_PROPERTIES",1))) {
    std::string final=parms.value_or_default("NUMBER_OF_WANG_LANDAU_STEPS","16");
    alps::RealVectorObsevaluator gf_=sim.get_measurements()["Coefficients "+final];
    std::valarray<double> gf=gf_.mean();
    alps::RealVectorObsevaluator m2_=sim.get_measurements()["Uniform Structure Factor Coefficients"];
    std::valarray<double> m2=m2_.mean();
    std::valarray<double> sm2;
    if (lattice.is_bipartite()) {
      alps::RealVectorObsevaluator sm2_=sim.get_measurements()["Staggered Structure Factor Coefficients"];
      std::valarray<double> sm2__=sm2_.mean();
      sm2.resize(sm2__.size());
      sm2=sm2__;
    }     
    alps::plot::Set<double> uniform_structure_factor_set;
    alps::plot::Set<double> uniform_susceptibility_set;
    alps::plot::Set<double> staggered_structure_factor_set;
    for (double T=T_min;T<T_max+T_delta/2.;T+=T_delta) {    
      double logbeta=-log(T);
      double maxx=gf[0];
      for (int i=1;i<gf.size();++i) // avoid overflow
        maxx=std::max(gf[i]+logbeta*i,maxx);
      double Z=0;
      double Sm2=0;
      double Ssm2=0;
      for (int i=0;i<gf.size();++i) {
        double factor=exp(gf[i]+logbeta*i-maxx);
        Z+=factor;
        Sm2+=m2[i]*factor;
        if (lattice.is_bipartite())
          Ssm2+=sm2[i]*factor;
      }
      double uniform_structure_factor=Sm2/Z/Ns;
      double uniform_susceptibility=uniform_structure_factor/T;
      uniform_structure_factor_set << T << uniform_structure_factor;
      uniform_susceptibility_set << T << uniform_susceptibility;
      if (lattice.is_bipartite()) {
        double staggered_structure_factor=Ssm2/Z/Ns;
        staggered_structure_factor_set << T <<  staggered_structure_factor;
      }
    }
    if (1) {
      alps::plot::Plot<double> my_plot("Uniform Structure Factor per Site versus Temperature",parms);
      my_plot.set_labels("Temperature","Uniform Structure Factor per Site");
      my_plot << uniform_structure_factor_set;
      alps::oxstream my_ox(ss+".plot.uniform_structure_factor.xml");
      my_ox << my_plot;
    }
    if (1) {
      alps::plot::Plot<double> my_plot("Uniform Susceptibility per Site versus Temperature",parms);
      my_plot.set_labels("Temperature","Uniform Susceptibility per Site");
      my_plot << uniform_susceptibility_set;
      alps::oxstream my_ox(ss+".plot.uniform_susceptibility.xml");
      my_ox << my_plot;
    }   
    if (lattice.is_bipartite()) {
      alps::plot::Plot<double> my_plot("Staggered Structure Factor per Site versus Temperature",parms);
      my_plot.set_labels("Temperature","Staggered Structure Factor per Site");
      my_plot << staggered_structure_factor_set;
      alps::oxstream my_ox(ss+".plot.staggered_structure_factor.xml");
      my_ox << my_plot;
    }
  }
}

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif
 
  alps::scheduler::SimpleMCFactory<alps::scheduler::DummyMCRun> factory;
  alps::scheduler::init(factory);

  int i=1;
  alps::Parameters parms;

  while (i<argc-1 && argv[i][0]=='-' && argv[i][1]=='-') {
    std::string name = argv[i]+2;
    if (name=="help") {
      std::cerr << "Usage: \n" << argv[0] << " [--T_MIN ...] [--T_MAX ...] [--DELTA_T ...] inputfile1 [inputfile2 [.....]] \n";
      ++i;
    }
    else {
      parms[name]=argv[i+1];
      i+=2;
    }
  }
  
  while (i < argc) {
    boost::filesystem::path p(argv[i]);
    evaluate(boost::filesystem::absolute(p),parms);
    ++i;
  }

 
#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  std::exit(-5);
}
#endif
  return 0;
}
