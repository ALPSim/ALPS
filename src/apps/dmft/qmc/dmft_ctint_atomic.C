// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// With ATOMIC=1 the interaction-expansion solver starts from the atomic-limit
// G0 and needs no INFILE, but its evaluation read G0 from INFILE (and MU for
// the high-frequency tails), so the run failed after simulating (#285). The
// evaluation now uses the same atomic G0. This test runs the Hubbard atom at
// half filling without INFILE and MU, in both measurement modes, and checks
// that G(tau) is the one of a run that reads G0 = 1/(i w_n) from INFILE with
// MU = 0, and the exact result. It also checks that INFILE, MU, H, EPS_f and
// EPSSQ_f given with ATOMIC=1 change nothing.

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <alps/solvers.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <string>
#include <vector>

namespace {
const double beta = 4.;
const double U = 2.;
const int n_tau = 50;
const int n_matsubara = 40;

// The solver appends to matrix_size and writes staggered_sz in the working
// directory.
struct working_directory {
  explicit working_directory(std::filesystem::path const & dir)
    : previous(std::filesystem::current_path()) { std::filesystem::current_path(dir); }
  ~working_directory() { std::filesystem::current_path(previous); }
  std::filesystem::path previous;
};

// G(tau) of one flavor of the Hubbard atom U (n_up - 1/2)(n_down - 1/2)
double exact_green(double tau) {
  return -0.5 * (std::exp(-tau * U / 2) + std::exp(-(beta - tau) * U / 2)) / (1 + std::exp(-beta * U / 2));
}

alps::params atomic_parameters(bool histogram) {
  alps::params p;
  p.supply("ATOMIC", "1");
  p.supply("BETA", std::to_string(beta));
  p.supply("U", std::to_string(U));
  p.supply("ALPHA", "-0.01");
  p.supply("N_TAU", std::to_string(n_tau));
  p.supply("N_MATSUBARA", std::to_string(n_matsubara));
  p.supply("HISTOGRAM_MEASUREMENT", histogram ? "1" : "0");
  p.supply("SWEEPS", "1000000");
  p.supply("THERMALIZATION", "100");
  p.supply("MAX_TIME", "60");
  p.supply("SEED", "5");
  return p;
}

// G0_f(i w_n) = 1/(i w_n), the atomic limit, or 1/(i w_n + 0.4)
std::string write_g0(std::filesystem::path const & dir, double shift) {
  const std::string file = (dir / ("g0_" + std::to_string(shift) + ".h5")).string();
  alps::hdf5::archive ar(file, "w");
  std::vector<std::complex<double> > g(n_matsubara);
  for (int n = 0; n < n_matsubara; ++n)
    g[n] = 1. / std::complex<double>(shift, (2 * n + 1) * std::acos(-1.) / beta);
  ar << alps::make_pvp("/G0_0", g);
  ar << alps::make_pvp("/G0_1", g);
  return file;
}

std::vector<double> green_tau(std::string const & file) {
  alps::hdf5::archive ar(file, "r");
  std::vector<double> g;
  ar >> alps::make_pvp("/G_tau/0/mean/value", g);
  return g;
}
}

TEST(DmftRegression, ctint_atomic_without_infile) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  for (bool histogram : {false, true}) {
    SCOPED_TRACE(histogram ? "HISTOGRAM_MEASUREMENT=1" : "HISTOGRAM_MEASUREMENT=0");
    const std::string output = (temporary.path() / (histogram ? "itime.h5" : "matsubara.h5")).string();
    alps::solvers::ctint(atomic_parameters(histogram), output);
    const std::vector<double> g = green_tau(output);
    ASSERT_EQ(g.size(), static_cast<std::size_t>(n_tau + 1));

    // the same G0 from INFILE, with the same random numbers
    alps::params p = atomic_parameters(histogram);
    p.supply("ATOMIC", "0");
    p.supply("INFILE", write_g0(temporary.path(), 0.));
    p.supply("MU", "0");
    p.supply("EPSSQ_0", "0");
    p.supply("EPSSQ_1", "0");
    const std::string reference = (temporary.path() / (histogram ? "itime_ref.h5" : "matsubara_ref.h5")).string();
    alps::solvers::ctint(p, reference);
    const std::vector<double> r = green_tau(reference);
    ASSERT_EQ(r.size(), g.size());
    double diff = 0;
    for (std::size_t t = 0; t < g.size(); ++t)
      diff = std::max(diff, std::abs(g[t] - r[t]));
    EXPECT_LT(diff, 1e-4);   // G0(tau) from INFILE is a Fourier transform on N_MATSUBARA frequencies

    // The largest error over 40 seeds is 0.019 in either measurement mode.
    for (int t = 0; t <= n_tau; t += 5)
      EXPECT_NEAR(g[t], exact_green(beta * t / n_tau), 0.05) << "tau index " << t;
  }
}

TEST(DmftRegression, ctint_atomic_ignores_infile) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  // a G0 the atomic run must not use
  const std::string g0 = write_g0(temporary.path(), 0.4);
  const std::string without = (temporary.path() / "without.h5").string();
  const std::string with = (temporary.path() / "with.h5").string();
  alps::solvers::ctint(atomic_parameters(false), without);
  alps::params p = atomic_parameters(false);
  p.supply("INFILE", g0);
  p.supply("MU", "0.7");
  p.supply("H", "0.3");
  p.supply("EPS_0", "0.2");
  p.supply("EPS_1", "-0.2");
  p.supply("EPSSQ_0", "5");
  p.supply("EPSSQ_1", "5");
  alps::solvers::ctint(p, with);
  EXPECT_EQ(green_tau(with), green_tau(without));
}
