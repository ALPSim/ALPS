/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2006-2009 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/parameter/parameters.h>
#include <boost/throw_exception.hpp>
#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <string>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
  try {
#endif

#ifndef BOOST_MSVC
  setenv("DIR", "/home/alps", 1);
#else
  _putenv("DIR=/home/alps");
#endif

  alps::Parameters params(std::cin);
  std::cout << params;

  // replace test
  params["L"] = 3;
  std::cout << "L = " << params["L"] << std::endl;

  // erase test
  params.erase("a0");
  std::cout << params;

  // test of copy constructor and copying a parameter from itself
  alps::Parameters params2(params);
  params2["N"] = params2["L"];
  std::cout << "N = " << params2["N"] << std::endl;

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
  return 0;
}
