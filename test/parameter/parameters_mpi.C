/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: parameters.C 2853 2008-06-17 13:59:59Z wistaria $ */

#include <alps/parameter/parameters.h>
#include <boost/config.hpp>
#include <boost/mpi.hpp>
#include <iostream>
#include <stdexcept>
#include <stdlib.h>

int main(int argc, char* argv[])
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;

  if (world.size() >= 2) {
    if (world.rank() == 0) {
      alps::Parameters params(std::cin);
      world.send(1, 0, params);
    } else if (world.rank() == 1) {
      alps::Parameters params;
      world.recv(0, 0, params);
      std::cout << params;
    }
  }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e) {
  std::cerr << "Caught exception: " << e.what() << "\n";
  exit(-1);
}
catch (...) {
  std::cerr << "Caught unknown exception\n";
  exit(-2);
}
#endif
}
