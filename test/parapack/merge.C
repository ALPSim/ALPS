/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/parapack/measurement.h>
#include <iomanip>
#include <iostream>

int main() {
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  std::cout << std::setprecision(3);

  alps::RealObservable obs_1("obs_1");
  for (int i = 0; i < 100; ++i)
    obs_1 << (double)i;
  std::cout << obs_1.mean() << " +/- " << obs_1.error() << std::endl;

  alps::ObservableSet obs_2;
  for (int i = 0; i < 100; ++i) {
    alps::ObservableSet obs;
    obs << alps::RealObservable("obs",10000);
    obs.reset(true);
    for (int j = 0; j < 100; ++j)
      obs["obs"] << (double)i;
    obs_2 << obs;
  }
  std::cout << dynamic_cast<alps::RealObsevaluator&>(obs_2["obs"]).mean() << " +/- "
            << dynamic_cast<alps::RealObsevaluator&>(obs_2["obs"]).error() << std::endl;
    

  alps::ObservableSet obs_3;
  for (int i = 0; i < 100; ++i) {
    alps::ObservableSet obs;
    obs << alps::RealObservable("obs");
    obs.reset(true);
    for (int j = 0; j < 100; ++j)
      obs["obs"] << (double)i;
    alps::merge_random_clone(obs_3, obs);
  }
  std::cout << dynamic_cast<alps::RealObservable&>(obs_3["obs"]).mean() << " +/- "
            << dynamic_cast<alps::RealObservable&>(obs_3["obs"]).error() << std::endl;
  return 0;
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
}
