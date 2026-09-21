/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parameter/parameters.h>
#include <alps/hdf5/archive.hpp>
#include <boost/filesystem/operations.hpp>
#include <iostream>
#include <cstdlib>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::string file = "parameters.h5";
  alps::Parameters params(std::cin);
  std::cout << params;
  {
    alps::hdf5::archive h5(file, "a");
    h5["/parameters"] << params;
  }
  params.clear();
  {
    alps::hdf5::archive h5(file, "a");
    h5["/parameters"] >> params;
  }
  std::cout << params;
  boost::filesystem::remove(boost::filesystem::path(file));

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exp) {
  std::cerr << exp.what() << std::endl;
  std::abort();
}
#endif
  return 0;
}
