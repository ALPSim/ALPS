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

#include "rng_helper.h"

namespace alps {

rng_helper::rng_helper(const Parameters& p) {
#ifdef ALPS_ENABLE_OPENMP_WORKER
  int nr = max_threads();
  engines_.resize(nr);
  generators_.resize(nr);
  #pragma omp parallel
  {
    int r = thread_id();
    engines_[r].reset(rng_factory.create(p.value_or_default("RNG", "mt19937")));
    generators_[r].reset(new generator_type(*engines_[r], boost::uniform_real<>()));
  }
#else
  engines_.resize(1);
  generators_.resize(1);
  engines_[0].reset(rng_factory.create(p.value_or_default("RNG", "mt19937")));
  generators_[0].reset(new generator_type(*engines_[0], boost::uniform_real<>()));
#endif
  init(p);
}

void rng_helper::init(const Parameters& p) {
  seed = static_cast<uint32_t>(p["WORKER_SEED"]);
  disorder_seed = static_cast<uint32_t>(p["DISORDER_SEED"]);
  alps::pseudo_des inigen(seed);
  for (unsigned int r = 0; r < engines_.size(); ++r) {
    engines_[r]->seed(inigen);
  }
  Disorder::seed(disorder_seed);
}

void rng_helper::load(IDump& dp) {
  std::string state;
  dp >> seed >> disorder_seed;
  for (unsigned int r = 0; r < engines_.size(); ++r) {
    dp >> state;
    std::stringstream rngstream(state);
    engines_[r]->read_all(rngstream);
  }
  Disorder::seed(disorder_seed);
}

void rng_helper::save(ODump& dp) const {
  dp << seed << disorder_seed;
  for (unsigned int r = 0; r < engines_.size(); ++r) {
    std::ostringstream rngstream;
    engines_[r]->write_all(rngstream);
    dp << rngstream.str();
  }
}

} // end namespace alps
