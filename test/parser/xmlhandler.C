/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/parser/xmlhandler.h>
#include <alps/parser/xmlparser.h>
#include <cstdlib>
#include <iostream>
#include <string>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  double v0;
  alps::SimpleXMLHandler<double> handler0("VALUE0", v0);

  double v1;
  alps::SimpleXMLHandler<double> handler1("VALUE1", v1, "value");

  alps::CompositeXMLHandler handler("TEST");
  handler.add_handler(handler0);
  handler.add_handler(handler1);

  alps::XMLParser parser(handler);
    
  parser.parse(std::cin);

  std::cout << v0 << std::endl
            << v1 << std::endl;

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exp) {
  std::cerr << exp.what() << std::endl;
  std::abort();
}
#endif
  return 0;
}
