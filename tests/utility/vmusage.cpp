/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010-2012 by Haruhiko Matsuo <halm@looper.t.u-tokyo.ac.jp>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/utility/vmusage.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#include <iostream>

int main(int argc, char** argv) {
  int pid = (argc == 1) ? -1 : boost::lexical_cast<int>(argv[1]);
  BOOST_FOREACH(alps::vmusage_type::value_type v, alps::vmusage(pid)) {
    std::cerr << v.first << " = " << v.second << "\n";
  }
  return 0;
}

  
