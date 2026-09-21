/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include "ising.h"
#include <alps/osiris/comm.h>

int main(int argc, char** argv)
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  int res=0;
  IsingFactory factory;
  
  alps::scheduler::SingleScheduler* s= alps::scheduler::start_single(factory,argc,argv);
  
  if (s) {
    // I'm the master and should actually do something
    alps::Parameters parms;
    std::cin >> parms;
    s->create_task(parms);
    res = s->run();
    std::cout << "Results: " << dynamic_cast<alps::scheduler::MCSimulation*>(s->get_task())->get_measurements();
    s->destroy_task();
  }
  
  alps::scheduler::stop_single();
  
  return res;
  
#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  alps::comm_exit(true);
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
}
