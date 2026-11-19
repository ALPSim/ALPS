// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// INSULATING starts the DMFT loop from the atomic G0 = 1/(i w_n + MU -+ H).
// It was enabled whenever the key was defined, so INSULATING = 0 also started
// from the atom, and it silently took precedence over G0OMEGA_INPUT.

#include <gtest/gtest.h>

#include <utility> // Boost.Math 1.76 includes this inside a namespace.
#include <boost/math/constants/constants.hpp>
#include "hilberttransformer.h"

#include <alps/parameter.h>

#include <complex>
#include <stdexcept>

namespace {
const double beta = 10., mu = 0.3, h = 0.1;
const unsigned int n_matsubara = 20;

alps::Parameters base_parameters() {
  alps::Parameters p;
  p["BETA"] = beta;
  p["MU"] = mu;
  p["H"] = h;
  p["t"] = 1;
  p["N"] = 50;
  p["NMATSUBARA"] = n_matsubara;
  p["FLAVORS"] = 2;
  return p;
}
}

TEST(DmftRegression, hilberttransformer_insulating) {
  alps::Parameters p = base_parameters();
  const matsubara_green_function_t noninteracting = GeneralFSHilbertTransformer(p, false).initial_G0(p);

  p["INSULATING"] = 1;
  const matsubara_green_function_t atomic = GeneralFSHilbertTransformer(p, false).initial_G0(p);
  p["INSULATING"] = 0;
  const matsubara_green_function_t off = GeneralFSHilbertTransformer(p, false).initial_G0(p);

  for (unsigned int n = 0; n < n_matsubara; ++n) {
    const std::complex<double> iw(0., (2 * n + 1) * boost::math::constants::pi<double>() / beta);
    for (unsigned int f = 0; f < 2; ++f) {
      const std::complex<double> expected = 1. / (iw + mu + (f % 2 ? h : -h));
      EXPECT_NEAR(std::abs(atomic(n, f) - expected), 0., 1e-12) << "n " << n << ", flavor " << f;
      EXPECT_EQ(off(n, f), noninteracting(n, f)) << "n " << n << ", flavor " << f;
    }
  }
}

TEST(DmftRegression, hilberttransformer_insulating_with_g0_input) {
  alps::Parameters p = base_parameters();
  p["INSULATING"] = 1;
  p["G0OMEGA_INPUT"] = "G0omega_input";
  EXPECT_THROW(GeneralFSHilbertTransformer(p, false).initial_G0(p), std::invalid_argument);
}

// The imaginary-time loop (OMEGA_LOOP = 0) reads its initial G0 from G0TAU_INPUT.
TEST(DmftRegression, hilberttransformer_insulating_with_g0tau_input) {
  alps::Parameters p = base_parameters();
  p["INSULATING"] = 1;
  p["G0TAU_INPUT"] = "G0tau_input";
  EXPECT_THROW(SemicircleHilbertTransformer(p).initial_G0(p), std::invalid_argument);
}
