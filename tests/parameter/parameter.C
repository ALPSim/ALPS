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

#include <alps/parameter/parameter.h>
#include <boost/throw_exception.hpp>
#include <cstdlib>
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

  std::string str;
  alps::Parameter p;
  while (std::getline(std::cin, str)) {
#ifndef BOOST_NO_EXCEPTIONS
    try {
#endif
    p.parse(str);
    str = p.value().c_str();
    std::cout << p.key() << " = ";
    if (str.find(' ') != std::string::npos)
      std::cout << '"' << str << '"';
    else
      std::cout << str;
    std::cout << ";\n";
#ifndef BOOST_NO_EXCEPTIONS
    }
    catch (std::exception& e) {
      std::cout << "Caught exception: " << e.what() << "\n";
    }
#endif
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
  return 0;
}
