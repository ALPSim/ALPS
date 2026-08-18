/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "mc_worker.h"

namespace alps {
namespace parapack {

//
// dumb_worker
//

dumb_worker::~dumb_worker() {}

void dumb_worker::print_copyright(std::ostream& out) {
  out << "ALPS/parapack dumb worker\n";
}

void dumb_worker::init_observables(Parameters const&, ObservableSet&) {}

void dumb_worker::run(ObservableSet&) {}

void dumb_worker::load(IDump&) {}

void dumb_worker::save(ODump&) const {}

bool dumb_worker::is_thermalized() const { return true; }

double dumb_worker::progress() const { return 1; }

//
// mc_worker
//

mc_worker::mc_worker(Parameters const& params)
  : abstract_worker(), rng_helper(params) {
}

mc_worker::~mc_worker() {}

void mc_worker::load_worker(IDump& dump) {
  abstract_worker::load_worker(dump);
  rng_helper::load(dump);
}

void mc_worker::save_worker(ODump& dump) const {
  abstract_worker::save_worker(dump);
  rng_helper::save(dump);
}

} // end namespace parapack
} // end namespace alps
