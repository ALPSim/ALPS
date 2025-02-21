/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2005-2008 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/parapack/integer_range.h>
#include <cstdlib>
#include <iostream>
#include <string>

int main()
{
  std::string str;
  alps::Parameters params;
  while (std::getline(std::cin, str) && str.size())
    params << alps::Parameter(str);
  std::cout << "Parameters:\n" << params;

  std::cout << "Test for integer_range<int>:\n";
  while (std::getline(std::cin, str)) {
    if (str.size() == 0 || str[0] == '\0') break;
    std::cout << "parse " << str << ": ";
    try {
      alps::integer_range<int> r(str, params);
      std::cout << "result " << r << std::endl;
    }
    catch (std::exception& exp) {
      std::cout << exp.what() << std::endl;
    }
  }
  std::cout << "Test for integer_range<unsigned int>:\n";
  while (std::getline(std::cin, str)) {
    if (!std::cin || str[0] == '\0') break;
    std::cout << "parse " << str << ": ";
    try {
      alps::integer_range<unsigned int> r(str, params);
      std::cout << "result " << r << std::endl;
    }
    catch (std::exception& exp) {
      std::cout << exp.what() << std::endl;
    }
  }

  alps::integer_range<int> r(0, 5);
  std::cout << "initial: " << r << std::endl;
  std::cout << "7 is included? " << r.is_included(7) << std::endl;
  r = 3;
  std::cout << "3 is assigned: " << r << std::endl;
  r.include(8);
  std::cout << "8 is included: " << r << std::endl;
  std::cout << "7 is included? " << r.is_included(7) << std::endl;

  alps::integer_range<int> s("[3:]");
  std::cout << "initial: " << s << std::endl;
  std::cout << "multiplied by 3.5: " << 3.5 * s << std::endl;
  std::cout << "multiplied by 2000000000: " << s * 2000000000 << std::endl;
  s *= 0.1;
  std::cout << "multiplied by 0.1: " << s << std::endl;

  alps::integer_range<int> t(0, 5);
  alps::integer_range<int> u(-2, 3);
  alps::integer_range<int> v(-2, 10);
  alps::integer_range<int> w(7, 10);
  std::cout << "initial: " << t << std::endl;
  std::cout << "overlap with " << u << ": " << overlap(t, u) << std::endl;
  std::cout << "union with " << u << ": " << unify(t, u) << std::endl;
  std::cout << "overlap with " << v << ": " << overlap(t, v) << std::endl;
  std::cout << "union with " << v << ": " << unify(t, v) << std::endl;
  std::cout << "overlap with " << w << ": " << overlap(t, w) << std::endl;
  std::cout << "union with " << w << ": ";
  try { std::cout << unify(t, w) << std::endl; }
  catch (std::exception& exp) { std::cout << exp.what() << std::endl; }
  return 0;
}
