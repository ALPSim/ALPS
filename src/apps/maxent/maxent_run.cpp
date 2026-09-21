/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2010 by Sebastian  Fuchs <fuchs@comp-phys.org>
*                       Thomas Pruschke <pruschke@comp-phys.org>
*                       Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "maxent.hpp"
#include <alps/solvers.hpp>

namespace {
bool stop_requested(boost::posix_time::ptime const& end_time) {
  static alps::ngs::signal signal;
  return !signal.empty() || boost::posix_time::second_clock::local_time() > end_time;
}
}

void alps::solvers::maxent(alps::params const& parms, std::string const& output_file) {
  MaxEntSimulation simulation(parms, output_file);
  simulation.run(boost::bind(&stop_requested, boost::posix_time::second_clock::local_time()
      + boost::posix_time::seconds(static_cast<int>(parms["MAX_TIME"] | 60))));
}
