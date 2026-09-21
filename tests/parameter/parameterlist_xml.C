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

#include <alps/parameter/parameterlist_p.h>
#include <iostream>
#include <cstdlib>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::ParameterList parameterlist;
  alps::ParameterListXMLHandler handler(parameterlist);
  
  alps::XMLParser parser(handler);
  parser.parse(std::cin);
  
  std::cout << "[Output to std::cout]\n";
  std::cout << parameterlist;
  
  std::cout << "\n[Output to alps::oxstream]\n";
  alps::oxstream oxs;
  oxs << parameterlist;

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exp) {
  std::cerr << exp.what() << std::endl;
  std::abort();
}
#endif
  return 0;
}
