/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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
