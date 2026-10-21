// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The interaction-expansion solver set its `thermalized` flag only in the
// constructor, so with THERMALIZATION > 0 fraction_completed() stayed 0, the
// run went on until MAX_TIME, and every step was measured. This test runs the
// solver with and without thermalization and checks that both stop after
// SWEEPS and keep SWEEPS/N_MEAS measurements, the same number for Sign and
// the signed observables.

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <alps/solvers.hpp>

#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

namespace {
const int sweeps = 2000;
const int n_meas = 25;
const int max_time = 10;

// The solver appends to matrix_size and writes staggered_sz in the working
// directory.
struct working_directory {
  explicit working_directory(std::filesystem::path const & dir)
    : previous(std::filesystem::current_path()) { std::filesystem::current_path(dir); }
  ~working_directory() { std::filesystem::current_path(previous); }
  std::filesystem::path previous;
};

std::uint64_t result_count(alps::hdf5::archive & ar, std::string const & name) {
  std::uint64_t count = 0;
  ar >> alps::make_pvp("/simulation/results/" + name + "/count", count);
  return count;
}

// Runs the solver on G0_f(i w_n) = 1/(i w_n + 0.4) for both flavors, which
// has a sign problem, and checks the number of measurements and the run time.
void check_run(std::filesystem::path const & dir, int thermalization) {
  SCOPED_TRACE("THERMALIZATION=" + std::to_string(thermalization));
  const double beta = 5.;
  const int n_matsubara = 32;
  const std::string g0 = (dir / "g0.h5").string();
  {
    alps::hdf5::archive ar(g0, "w");
    std::vector<std::complex<double> > g(n_matsubara);
    for (int n = 0; n < n_matsubara; ++n)
      g[n] = 1. / std::complex<double>(0.4, (2 * n + 1) * std::acos(-1.) / beta);
    ar << alps::make_pvp("/G0_0", g);
    ar << alps::make_pvp("/G0_1", g);
  }

  alps::params p;
  p.supply("INFILE", g0);
  p.supply("BETA", std::to_string(beta));
  p.supply("U", "4");
  p.supply("MU", "1");
  p.supply("ALPHA", "0.3");
  p.supply("N_ORBITALS", "2");
  p.supply("N_MATSUBARA", std::to_string(n_matsubara));
  p.supply("N_TAU", "64");
  p.supply("HISTOGRAM_MEASUREMENT", "1");
  p.supply("SWEEPS", std::to_string(sweeps));
  p.supply("THERMALIZATION", std::to_string(thermalization));
  p.supply("N_MEAS", std::to_string(n_meas));
  p.supply("MAX_TIME", std::to_string(max_time));
  p.supply("SEED", "3");

  const std::string output = (dir / ("out" + std::to_string(thermalization) + ".h5")).string();
  const auto start = std::chrono::steady_clock::now();
  alps::solvers::ctint(p, output);
  const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
  EXPECT_LT(elapsed.count(), max_time / 2.) << "the run did not stop after SWEEPS";

  alps::hdf5::archive ar(output, "r");
  const std::uint64_t sign = result_count(ar, "Sign");
  EXPECT_EQ(sign, static_cast<std::uint64_t>(sweeps / n_meas));
  EXPECT_EQ(result_count(ar, "Sz_0"), sign);
  EXPECT_EQ(result_count(ar, "W_0_0_0"), sign);
  EXPECT_EQ(result_count(ar, "PertOrder"), sign);
}
}

TEST(DmftRegression, ctint_thermalization_ends) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  check_run(temporary.path(), 0);
  check_run(temporary.path(), 100);
}
