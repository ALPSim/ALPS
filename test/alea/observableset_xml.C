/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2006 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <iostream>
#include <alps/alea/observableset_p.h>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::ObservableSet obs;
  alps::ObservableSetXMLHandler handler(obs);
  alps::XMLParser parser(handler);
  parser.parse(std::cin);
  
  alps::oxstream oxs;
  obs.write_xml(oxs);

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exc) {
  std::cerr << exc.what() << "\n";
  return -1;
}
catch (...) {
  std::cerr << "Fatal Error: Unknown Exception!\n";
  return -2;
}
#endif
  return 0;
}
