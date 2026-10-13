/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/parameter/parameters_p.h>
#include <alps/osiris/xdrdump.h>
#include <iostream>
#include <cstdlib>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::Parameters parameters;
  alps::ParametersXMLHandler handler(parameters);
  
  alps::XMLParser parser(handler);
  parser.parse(std::cin);
  
  std::cout << parameters;

  {
    alps::OXDRFileDump od(boost::filesystem::path("parameters.dump"));
    od << parameters;
  }

  parameters.clear();
  
  {
    alps::IXDRFileDump id(boost::filesystem::path("parameters.dump"));
    id >> parameters;
  }

  std::cout << parameters;

  alps::oxstream oxs;
  oxs << parameters;

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exp) {
  std::cerr << exp.what() << std::endl;
  std::abort();
}
#endif
  return 0;
}
