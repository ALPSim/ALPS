/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2011 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/parapack/exchange.h>
#include <algorithm>
#include <iostream>
#include <boost/random.hpp>

int main() {
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  typedef alps::parapack::exmc::inverse_temperature_set inverse_temperature_set;

  boost::mt19937 eng;
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
    rng(eng, boost::uniform_real<>());

  std::cout << std::setprecision(3);

  int n = 10;

  std::vector<double> x, y;
  for (int i = 0; i < n; ++i) {
    x.push_back(rng());
    y.push_back(rng());
  }
  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  for (int i = 0; i < n; ++i) {
    std::cout << x[i] << '\t' << y[i] << std::endl;
  }


  double a, b;

  boost::tie(a,b) = inverse_temperature_set::linear_regression(10, 0, x, y);
  std::cout << "#\t" << 10 << '\t' << 0 << '\t' << a << '\t' << b << std::endl;

  boost::tie(a,b) = inverse_temperature_set::linear_regression(2, 2, x, y);
  std::cout << "#\t" << 2 << '\t' << 2 << '\t' << a << '\t' << b << std::endl;

  boost::tie(a,b) = inverse_temperature_set::linear_regression(5, 3, x, y);
  std::cout << "#\t" << 5 << '\t' << 3 << '\t' << a << '\t' << b << std::endl;

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
