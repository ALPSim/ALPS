/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2009 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/parameter/parameterlist.h>
#include <alps/osiris/xdrdump.h>
#include <boost/filesystem/operations.hpp>
#include <boost/throw_exception.hpp>
#include <stdlib.h>
#include <iostream>
#include <stdexcept>

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

  boost::filesystem::path path("parameterlist.dump");

  alps::ParameterList params(std::cin);
  std::cout << params;

  {
    alps::OXDRFileDump od(path);
    od << params;
  }

  params.clear();

  {
    alps::IXDRFileDump id(path);
    id >> params;
  }

  std::cout << params;
  boost::filesystem::remove(path);

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& e)
{
  std::cerr << "Caught exception: " << e.what() << "\n";
  exit(-1);
}
catch (...)
{
  std::cerr << "Caught unknown exception\n";
  exit(-2);
}
#endif
  return 0;
}
