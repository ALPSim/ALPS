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

#include <alps/parapack/temperature_scan.h>
#include <iomanip>
#include <iostream>

class my_worker {
public:
  my_worker(alps::Parameters const&) {}
  void init_observables(alps::Parameters const&, alps::ObservableSet const&) {}
  void run(alps::ObservableSet const&) {}
  void set_beta(double beta) { std::cout << "T = " << 1/beta; }
  void save(alps::ODump&) const {}
  void load(alps::IDump&) {}
};

int main() {
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif
  std::cout << std::setprecision(3);

  alps::Parameters params(std::cin);
  std::vector<alps::ObservableSet> obs;
  alps::parapack::temperature_scan_adaptor<my_worker> worker(params);
  worker.init_observables(params, obs);

  int count = 0;
  while (worker.progress() < 1) {
    std::cout << count++ << ", ";
    worker.run(obs);
    std::cout << ", progress = " << worker.progress() << std::endl;
  }

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
