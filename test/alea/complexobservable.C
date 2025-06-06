/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Synge Todo <wistaria@comp-phys.org>
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

/* $Id: detailedbinning.C 3654 2010-01-06 23:47:46Z troyer $ */

#include <iostream>
#include <alps/alea.h>
#include <boost/random.hpp>
#include <alps/alea/abstractsimpleobservable.ipp>
#include <alps/alea/simpleobservable.ipp>
#include <alps/alea/simpleobseval.ipp>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  int n = 1000;

  // random number generator
  typedef boost::minstd_rand0 random_base_type;
  typedef boost::uniform_01<random_base_type> random_type;
  random_base_type random_int;
  random_type random(random_int);

  // observable
  alps::ComplexObservable obs_a("Complex Observable A");
  alps::ComplexObservable obs_b("Complex Observable B");

  for(int i = 0; i < n; ++i) {
    double re = random();
    double im = random()+1;
    obs_a << std::complex<double>(re, im);
    re = random()+3;
    im = random()+1;
    obs_b << std::complex<double>(re, im);
  }

  // output
  std::cout << obs_a << obs_b;

  alps::ComplexObsevaluator eval_a = obs_a;
  alps::ComplexObsevaluator eval_b = obs_b;
  alps::ComplexObsevaluator prod = eval_a * eval_b;
  
  std::cout << prod;
  // std::cout << real(prod); // not yet implemented
  // std::cout << imag(prod); // not yet implemented
  // std::cout << abs(prod);  // not yet implemented
  // std::cout << arg(prod);  // not yet implemented

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
