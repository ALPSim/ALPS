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
#include <alps/utilities/signal.hpp>
#include <alps/utilities/copyright.hpp>
#include <alps/cli.hpp>
#include <alps/solvers.hpp>
#include <boost/bind/bind.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace {
bool stop_requested(boost::posix_time::ptime const& end_time) {
  static alps::signal signal;
  return !signal.empty() || boost::posix_time::second_clock::local_time() > end_time;
}
}

void alps::solvers::maxent(alps::params const& input, std::string const& output_file) {
  alps::params parms(input);
  MaxEntSimulation::define_parameters(parms);
  std::ostringstream missing;
  if (parms.has_missing(missing))
    throw std::invalid_argument("maxent: " + missing.str());
  MaxEntSimulation simulation(parms, output_file);
  if (alps::cli_is_master()) {
    alps::print_copyright(std::cout, "maxent");
    std::cout << std::flush;
    parms.has_unused(std::cerr);
  }
  simulation.run(boost::bind(&stop_requested, boost::posix_time::second_clock::local_time()
      + boost::posix_time::seconds(static_cast<int>(parms["MAX_TIME"]))));
}
