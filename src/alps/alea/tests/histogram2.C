/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2013 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// test for merging two empty histogram observables

#include <alps/alea.h>

int main() {
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  alps::oxstream oxs;

  alps::ObservableSet obs1;
  obs1 << alps::HistogramObservable<alps::int32_t>("histogram 1", 0, 10);
  obs1.write_xml(oxs);

  alps::ObservableSet obs2;
  obs2 << alps::HistogramObservable<alps::int32_t>("histogram 1", 0, 10);
  obs2.write_xml(oxs);

  obs1 << obs2;
  obs1.write_xml(oxs);

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
