/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2006-2009 by Synge Todo <wistaria@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
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
