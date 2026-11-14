// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The CT-INT solver that dmft runs in-process (SOLVER = "Interaction
// Expansion") checked only whether ATOMIC was defined, so ATOMIC = 0 also
// replaced the bath by the atomic limit (#347). This test calls the solver as
// dmft does, with ATOMIC = 0 and the bath G0 = 1/(i w_n + MU) of a Hubbard atom
// away from half filling, and compares G(tau) with the exact result.

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include "interaction_expansion/interaction_expansion.hpp"

#include <alps/parameter.h>

#include <cmath>
#include <complex>
#include <filesystem>

namespace {
const double beta = 4.;
const double U = 2.;
const double mu = 0.4;
const unsigned int n_tau = 50;
const unsigned int n_matsubara = 40;

// The solver writes BASENAME.h5 and matrix_size in the working directory.
struct working_directory {
  explicit working_directory(std::filesystem::path const & dir)
    : previous(std::filesystem::current_path()) { std::filesystem::current_path(dir); }
  ~working_directory() { std::filesystem::current_path(previous); }
  std::filesystem::path previous;
};

// G(tau) of one flavor of the atom U (n_up - 1/2)(n_down - 1/2) - MU (n_up + n_down)
double exact_green(double tau) {
  const double e0 = U / 4, e1 = -U / 4 - mu, e2 = U / 4 - 2 * mu;
  const double z = std::exp(-beta * e0) + 2 * std::exp(-beta * e1) + std::exp(-beta * e2);
  return -(std::exp(-beta * e0 + tau * (e0 - e1)) + std::exp(-beta * e1 + tau * (e1 - e2))) / z;
}
}

TEST(DmftRegression, interaction_expansion_atomic_zero) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());

  alps::Parameters p;
  p["BASENAME"] = "atomic_zero";
  p["FLAVORS"] = 2;
  p["SITES"] = 1;
  p["BETA"] = beta;
  p["U"] = U;
  p["MU"] = mu;
  p["H"] = 0;
  p["EPSSQ_0"] = 0;
  p["EPSSQ_1"] = 0;
  p["ALPHA"] = -0.01;
  p["N"] = n_tau;
  p["NMATSUBARA"] = n_matsubara;
  p["SWEEPS"] = 1000000;
  p["THERMALIZATION"] = 100;
  p["MAX_TIME"] = 60;
  p["SEED"] = 5;
  p["ATOMIC"] = 0;

  matsubara_green_function_t bath(n_matsubara, 1, 2);
  for (unsigned int f = 0; f < 2; ++f)
    for (unsigned int n = 0; n < n_matsubara; ++n)
      bath(n, 0, 0, f) = 1. / std::complex<double>(mu, (2 * n + 1) * M_PI / beta);

  char name[] = "dmft_interaction_expansion_atomic";
  char *argv[] = {name, nullptr};
  alps::scheduler::BasicFactory<InteractionExpansionSim, HubbardInteractionExpansionRun> factory;
  alps::ImpuritySolver solver(factory, p, 1, argv);
  const itime_green_function_t g = solver.solve_omega(bath, p).second;

  // Over 12 seeds the largest error is 0.018; with the atomic limit it was 0.05.
  for (unsigned int f = 0; f < 2; ++f)
    for (unsigned int t = 0; t <= n_tau; ++t)
      EXPECT_NEAR(g(t, 0, 0, f), exact_green(beta * t / n_tau), 0.03)
        << "flavor " << f << ", tau index " << t;
}
