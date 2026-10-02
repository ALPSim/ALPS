/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/clone_timer.h>
#include <iostream>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  double progress = 0;
  alps::clone_timer timer(boost::posix_time::milliseconds(100), progress);

  alps::clone_timer::loops_t loops = 1024;
  alps::clone_timer::time_t start = timer.current_time();
  alps::clone_timer::time_t end = start + boost::posix_time::seconds(1);
  std::cerr << "start time = " << start << std::endl;

  while (true) {
    std::cerr << "loops = " << loops << std::endl;
    for (alps::clone_timer::loops_t i = 0; i < loops; ++i) {
      timer.current_time();
    }
    alps::clone_timer::time_t current = timer.current_time();
    std::cerr << "current time = " << current << " (" << (current - start) << ")\n";
    if (current > end) break;
    loops = timer.next_loops(loops, current);
  }
  end += boost::posix_time::seconds(1);
  while (true) {
    std::cerr << "loops = " << loops << std::endl;
    for (alps::clone_timer::loops_t i = 0; i < loops; ++i) {
      timer.current_time();
      timer.current_time();
      timer.current_time();
      timer.current_time();
    }
    alps::clone_timer::time_t current = timer.current_time();
    std::cerr << "current time = " << current << " (" << (current - start) << ")\n";
    if (current > end) break;
    loops = timer.next_loops(loops, current);
  }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exp) {
  std::cerr << exp.what() << std::endl;
  std::abort();
}
#endif
  return 0;
}
