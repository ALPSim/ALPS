/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/parapack/exp_number.h>
#include <alps/parser/xmlstream.h>
#include <iostream>

int main() {
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::exp_double x = 3;
  alps::exp_double y = 2;
  alps::exp_double z = x + y;
  alps::exp_double v = alps::exp_value(10000);
  alps::exp_double w = alps::exp_value(9999);
  alps::exp_double p = -2.5;

  std::cout << "x = " << x << " = " << alps::precision(static_cast<double>(x), 6) << std::endl;
  std::cout << "y = " << y << " = " << alps::precision(static_cast<double>(y), 6) << std::endl;
  std::cout << "x + y = " << x + y << " = " << alps::precision(static_cast<double>(x + y), 6) << std::endl;
  // std::cout << "x - y = " << x - y << " = " << static_cast<double>(x - y) << std::endl;
  std::cout << "x * y = " << x * y << " = " << alps::precision(static_cast<double>(x * y), 6) << std::endl;
  std::cout << "x / y = " << x / y << " = " << alps::precision(static_cast<double>(x / y), 6) << std::endl;
  std::cout << "x + 1.2 = " << x + 1.2 << " = " << alps::precision(static_cast<double>(x + 1.2), 6) << std::endl;
  std::cout << "x - 1.2 = " << x - 1.2 << " = " << alps::precision(static_cast<double>(x - 1.2), 6) << std::endl;
  std::cout << "x - 3.5 = " << x - 3.5 << " = " << alps::precision(static_cast<double>(x - 3.5), 6) << std::endl;
  std::cout << "x * 1.2 = " << x * 1.2 << " = " << alps::precision(static_cast<double>(x * 1.2), 6) << std::endl;
  std::cout << "x * 2 = " << x * 2 << " = " << alps::precision(static_cast<double>(x * 2), 6) << std::endl;
  std::cout << "x / 1.2 = " << x / 1.2 << " = " << alps::precision(static_cast<double>(x / 1.2), 6) << std::endl;
  std::cout << "x / 3 = " << x / 3 << " = " << alps::precision(static_cast<double>(x / 3), 6) << std::endl;
  std::cout << "3.5 + x = " << 3.5 + x << " = " << alps::precision(static_cast<double>(3.5 + x), 6) << std::endl;
  std::cout << "3 + x = " << 3 + x << " = " << alps::precision(static_cast<double>(3 + x), 6) << std::endl;
  std::cout << "3.5 - x = " << 3.5 - x << " = " << alps::precision(static_cast<double>(3.5 - x), 6) << std::endl;
  std::cout << "5 - x = " << 5 - x << " = " << alps::precision(static_cast<double>(4 - x), 6) << std::endl;
  std::cout << "1.2 - x = " << 1.2 - x << " = " << alps::precision(static_cast<double>(1.2 - x), 6) << std::endl;
  std::cout << "3.5 * x = " << 3.5 * x << " = " << alps::precision(static_cast<double>(3.5 * x), 6) << std::endl;
  std::cout << "3 * x = " << 3 * x << " = " << alps::precision(static_cast<double>(3 * x), 6) << std::endl;
  std::cout << "3.5 / x = " << 3.5 / x << " = " << alps::precision(static_cast<double>(3.5 / x), 6) << std::endl;
  std::cout << "3 / x = " << 3 / x << " = " << alps::precision(static_cast<double>(3 / x), 6) << std::endl;

  std::cout << "v = " << v << " = " << alps::precision(static_cast<double>(v), 6) << std::endl;
  std::cout << "w = " << w << " = " << alps::precision(static_cast<double>(w), 6) << std::endl;
  std::cout << "pow(v,3) = " << pow(v,3) << " = " << alps::precision(static_cast<double>(pow(v,3)), 6)
            << std::endl;
  std::cout << "v + w = " << v + w << " = " << alps::precision(static_cast<double>(v + w), 6) << std::endl;
  std::cout << "v - w = " << v - w << " = " << alps::precision(static_cast<double>(v - w), 6) << std::endl;
  std::cout << "v * w = " << v * w << " = " << alps::precision(static_cast<double>(v * w), 6) << std::endl;
  std::cout << "v / w = " << v / w << " = " << alps::precision(static_cast<double>(v / w), 6) << std::endl;
  std::cout << "v > w = " << (v > w) << std::endl;
  std::cout << "v >= w = " << (v >= w) << std::endl;
  std::cout << "v == w = " << (v == w) << std::endl;
  std::cout << "v <= w = " << (v <= w) << std::endl;
  std::cout << "v < w = " << (v < w) << std::endl;
  std::cout << "v > v = " << (v > v) << std::endl;
  std::cout << "v >= v = " << (v >= v) << std::endl;
  std::cout << "v == v = " << (v == v) << std::endl;
  std::cout << "v <= v = " << (v <= v) << std::endl;
  std::cout << "v < v = " << (v < v) << std::endl;

  std::cout << "p = " << p << " = " << alps::precision(static_cast<double>(p), 6) << std::endl;
  std::cout << "-p = " << -p << " = " << alps::precision(static_cast<double>(-p), 6) << std::endl;
  std::cout << "p - p = " << p - p << " = " << alps::precision(static_cast<double>(p - p), 6) << std::endl;
  std::cout << "x + p = " << x + p << " = " << alps::precision(static_cast<double>(x + p), 6) << std::endl;
  std::cout << "x - p = " << x - p << " = " << alps::precision(static_cast<double>(x - p), 6) << std::endl;
  std::cout << "x * p = " << x * p << " = " << alps::precision(static_cast<double>(x * p), 6) << std::endl;
  std::cout << "x / p = " << x / p << " = " << alps::precision(static_cast<double>(x / p), 6) << std::endl;
  std::cout << "y + p = " << y + p << " = " << alps::precision(static_cast<double>(y + p), 6) << std::endl;
  std::cout << "y - p = " << y - p << " = " << alps::precision(static_cast<double>(y - p), 6) << std::endl;
  std::cout << "y * p = " << y * p << " = " << alps::precision(static_cast<double>(y * p), 6) << std::endl;
  std::cout << "y / p = " << y / p << " = " << alps::precision(static_cast<double>(y / p), 6) << std::endl;
  std::cout << "x > p = " << (x > p) << std::endl;
  std::cout << "x >= p = " << (x >= p) << std::endl;
  std::cout << "x == p = " << (x == p) << std::endl;
  std::cout << "x <= p = " << (x <= p) << std::endl;
  std::cout << "x < p = " << (x < p) << std::endl;
  std::cout << "y > p = " << (y > p) << std::endl;
  std::cout << "y >= p = " << (y >= p) << std::endl;
  std::cout << "y == p = " << (y == p) << std::endl;
  std::cout << "y <= p = " << (y <= p) << std::endl;
  std::cout << "y < p = " << (y < p) << std::endl;
  std::cout << "p == p = " << (p == p) << std::endl;
  std::cout << "p == -p = " << (p == -p) << std::endl;
  std::cout << "x > 2 = " << (x > 2) << std::endl;
  std::cout << "x >= 2 = " << (x >= 2) << std::endl;
  std::cout << "x == 2 = " << (x == 2) << std::endl;
  std::cout << "x <= 2 = " << (x <= 2) << std::endl;
  std::cout << "x < 2 = " << (x < 2) << std::endl;

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
