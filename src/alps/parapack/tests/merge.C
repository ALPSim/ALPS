/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
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
