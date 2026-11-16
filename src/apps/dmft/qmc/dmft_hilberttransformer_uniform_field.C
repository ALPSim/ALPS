// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// A uniform field H shifts flavor 0 by MU - H and flavor 1 by MU + H, as the
// impurity solvers do. Without ANTIFERROMAGNET, the general Hilbert transform
// gave both flavors MU - H, and the semicircle transform coupled each flavor
// to the other flavor's G, which is the antiferromagnetic self-consistency.

#include <gtest/gtest.h>

#include <utility> // Boost.Math 1.76 includes this inside a namespace.
#include <boost/math/constants/constants.hpp>
#include "hilberttransformer.h"

#include <alps/parameter.h>

#include <complex>

namespace {
const double beta = 10., mu = 0.3, h = 0.2, t = 1.;
const unsigned int n_matsubara = 20;

std::complex<double> matsubara(unsigned int n) {
  return std::complex<double>(0., (2 * n + 1) * boost::math::constants::pi<double>() / beta);
}

// Local G of the semicircular DOS with second moment t^2 at complex energy z.
std::complex<double> semicircle_green(std::complex<double> z) {
  return z * (1. - std::sqrt(1. - 4. * t * t / (z * z))) / (2. * t * t);
}

double field(unsigned int flavor) { return flavor % 2 ? h : -h; }

alps::Parameters parameters() {
  alps::Parameters p;
  p["BETA"] = beta;
  p["t"] = t;
  p["FLAVORS"] = 2;
  p["ANTIFERROMAGNET"] = "false";
  return p;
}
}

// With no self-energy the transform returns the non-interacting local G.
TEST(DmftRegression, hilberttransformer_uniform_field_general) {
  alps::Parameters p = parameters();
  GeneralFSHilbertTransformer transform(p, false);

  matsubara_green_function_t G(n_matsubara, 1, 2), G0(n_matsubara, 1, 2);
  for (unsigned int n = 0; n < n_matsubara; ++n)
    for (unsigned int f = 0; f < 2; ++f) G(n, f) = G0(n, f) = 1. / (matsubara(n) + mu);

  const matsubara_green_function_t out = transform(G, G0, mu, h, beta);
  for (unsigned int n = 0; n < n_matsubara; ++n)
    for (unsigned int f = 0; f < 2; ++f)
      EXPECT_NEAR(std::abs(out(n, f) - semicircle_green(matsubara(n) + mu + field(f))), 0., 1e-12)
        << "n " << n << ", flavor " << f;
}

// The Bethe-lattice condition G0_f^-1 = i w_n + MU -+ H - t^2 G_f, flavor by flavor.
TEST(DmftRegression, hilberttransformer_uniform_field_semicircle) {
  alps::Parameters p = parameters();
  SemicircleFSHilbertTransformer transform(p);

  matsubara_green_function_t G(n_matsubara, 1, 2), unused(n_matsubara, 1, 2);
  for (unsigned int n = 0; n < n_matsubara; ++n) {
    G(n, 0) = 1. / (matsubara(n) + 0.5);
    G(n, 1) = 1. / (matsubara(n) - 0.7);
  }

  const matsubara_green_function_t out = transform(G, unused, mu, h, beta);
  for (unsigned int n = 0; n < n_matsubara; ++n)
    for (unsigned int f = 0; f < 2; ++f)
      EXPECT_NEAR(std::abs(out(n, f) - 1. / (matsubara(n) + mu + field(f) - t * t * G(n, f))), 0., 1e-12)
        << "n " << n << ", flavor " << f;
}
