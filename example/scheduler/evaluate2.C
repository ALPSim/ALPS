/***************************************************************************
* ALPS++/scheduler library
*
* scheduler/convert2xml.C   convert old scheduler files to XML
*
* $Id$
*
* Copyright (C) 2002 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
**************************************************************************/

#include <alps/scheduler.h>
#include <fstream>

void evaluate(const boost::filesystem::path& p, std::ostream& out) {
  alps::ProcessList nowhere;
  alps::scheduler::MCSimulation sim(nowhere,p);
  alps::RealObsevaluator m2=sim.get_measurements()["Magnetization^2"];
  alps::RealObsevaluator m4=sim.get_measurements()["Magnetization^4"];
#ifdef ALPS_HAVE_VALARRAY
  alps::RealVectorObsevaluator corr=sim.get_measurements()["Correlations"];
#endif
  alps::RealObsevaluator binder=1.-3.*m4/(m2*m2);
 
  binder.rename("Real binder cumulant of Magnetization");
#ifdef ALPS_HAVE_VALARRAY
  out << corr << "\n";
#endif
  out << binder << "\n";
  binder.write_xml(out);
  sim << binder;
  sim.checkpoint(p);
}

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::scheduler::SimpleFactory<alps::scheduler::DummyMCRun> factory;
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
