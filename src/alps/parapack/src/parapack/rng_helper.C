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
#include <alps/parapack/rng_helper.h>
#include <alps/input/parapack.hpp>

namespace alps {

rng_helper::rng_helper(const Parameters& p) {
  input::parapack_rng_parameters const params = input::read_parapack_rng_parameters(p);
  std::string const& rng = params.rng;
#ifdef ALPS_ENABLE_OPENMP
  int nr = max_threads();
  engines_.resize(nr);
  generators_.resize(nr);
  #pragma omp parallel
  {
    int r = thread_id();
    engines_[r].reset(rng_factory.create(rng));
    generators_[r].reset(new generator_type(*engines_[r], boost::uniform_real<>()));
  }
#else
  engines_.resize(1);
  generators_.resize(1);
  engines_[0].reset(rng_factory.create(rng));
  generators_[0].reset(new generator_type(*engines_[0], boost::uniform_real<>()));
#endif
  init(params);
}

void rng_helper::init(input::parapack_rng_parameters const& p) {
  seed = p.worker_seed;
  disorder_seed = p.disorder_seed;
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
