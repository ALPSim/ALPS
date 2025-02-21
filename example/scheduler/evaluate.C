/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>
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
