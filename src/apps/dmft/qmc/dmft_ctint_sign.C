// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The interaction-expansion solver with a sign problem, against the exact
// Hubbard atom U (n_up - 1/2)(n_down - 1/2) at half filling (#390):
//
// - Every observable is accumulated as value times sign. Only Sz_i, Sz2_i and
//   Sz0_Szi were divided by the mean sign, so G and the densities came out
//   multiplied by <sign>.
// - reset_perturbation_series set the sign to 1 every RECALC_PERIOD steps
//   without changing the configuration, which corrupted every signed average.
//   G(tau) of the atom is the same for U and -U and so cannot see this; the
//   double occupancy and <Sz^2> can.
// - Sz2_i was (n_up - n_down)^2 of the configuration's densities, which leaves
//   out the n_f (1 - n_f) of each flavor that Wick's theorem gives for n_f^2.
//
// U = -2 with ALPHA = 0.3 has <sign> ~ 0.5; U = 2 with ALPHA = -0.01 has none.

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <alps/solvers.hpp>

#include <cmath>
#include <filesystem>
#include <string>
#include <vector>

#ifdef ALPS_HAVE_MPI
#include <mpi.h>
#endif

namespace {
// Under MPI only rank 0 writes the result file.
bool is_master() {
#ifdef ALPS_HAVE_MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank == 0;
#else
  return true;
#endif
}

const double beta = 4.;
const int n_tau = 50;
// A deviation of more than this many error bars fails. Over seeds 5, 11, 23,
// 37 and 41 the largest deviation is 2.0 error bars.
const double n_sigma = 4.;

struct working_directory {
  explicit working_directory(std::filesystem::path const & dir)
    : previous(std::filesystem::current_path()) { std::filesystem::current_path(dir); }
  ~working_directory() { std::filesystem::current_path(previous); }
  std::filesystem::path previous;
};

// G(tau) of one flavor; the same for U and -U
double exact_green(double U, double tau) {
  const double a = std::abs(U) / 2;
  return -0.5 * (std::exp(-tau * a) + std::exp(-(beta - tau) * a)) / (1 + std::exp(-beta * a));
}
// <(n_up - n_down)^2>, the weight of the singly occupied states
double exact_sz2(double U) { return 1. / (1. + std::exp(-beta * U / 2)); }
// <n_up n_down> = (1 - <Sz^2>) / 2 at half filling
double exact_double_occupancy(double U) { return (1. - exact_sz2(U)) / 2; }

alps::params parameters(double U, double alpha, bool histogram) {
  alps::params p;
  p.supply("ATOMIC", "1");
  p.supply("BETA", std::to_string(beta));
  p.supply("U", std::to_string(U));
  p.supply("ALPHA", std::to_string(alpha));
  p.supply("N_TAU", std::to_string(n_tau));
  p.supply("N_MATSUBARA", "40");
  p.supply("HISTOGRAM_MEASUREMENT", histogram ? "1" : "0");
  p.supply("SWEEPS", "20000000");
  p.supply("THERMALIZATION", "1000");
  p.supply("MAX_TIME", "120");
  p.supply("SEED", "5");
  return p;
}

struct estimate { double mean, error; };

template <class T> T read(std::string const & file, std::string const & path) {
  alps::hdf5::archive ar(file, "r");
  T value;
  ar >> alps::make_pvp(path, value);
  return value;
}

estimate scalar_result(std::string const & file, std::string const & name) {
  const std::string path = "/simulation/results/" + name + "/mean/";
  return {read<double>(file, path + "value"), read<double>(file, path + "error")};
}

void expect_exact(estimate const & x, double exact, char const * what) {
  EXPECT_GT(x.error, 0.) << what;
  EXPECT_LT(x.error, 0.01) << what;
  EXPECT_LT(std::abs(x.mean - exact), n_sigma * x.error)
      << what << " = " << x.mean << " +- " << x.error << ", exact " << exact;
}
}

// Matsubara measurement: G, the densities and the double occupancy
TEST(DmftRegression, ctint_sign_matsubara) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  const double U = -2.;
  const std::string output = (temporary.path() / "matsubara.h5").string();
  alps::solvers::ctint(parameters(U, 0.3, false), output);
  if (!is_master())
    return;

  const double sign = scalar_result(output, "Sign").mean;
  EXPECT_LT(sign, 0.8);   // the test needs a sign problem
  EXPECT_GT(sign, 0.2);

  const std::string dens = "/simulation/results/densities/mean/";
  const std::vector<double> n = read<std::vector<double> >(output, dens + "value");
  const std::vector<double> n_err = read<std::vector<double> >(output, dens + "error");
  ASSERT_EQ(n.size(), 2u);
  for (int z = 0; z < 2; ++z)
    expect_exact({n[z], n_err[z]}, 0.5, "density");
  expect_exact(scalar_result(output, "density_correlation"), exact_double_occupancy(U),
               "double occupancy");

  // G carries no error bars.
  const std::vector<double> g = read<std::vector<double> >(output, "/G_tau/0/mean/value");
  ASSERT_EQ(g.size(), static_cast<std::size_t>(n_tau + 1));
  for (int t = 0; t <= n_tau; t += 5)
    EXPECT_NEAR(g[t], exact_green(U, beta * t / n_tau), 0.01) << "tau index " << t;
}

// Imaginary-time measurement: <Sz^2> with and without a sign problem
TEST(DmftRegression, ctint_sign_sz2) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  struct run { double U, alpha; };
  for (run r : {run{-2., 0.3}, run{2., -0.01}}) {
    SCOPED_TRACE("U = " + std::to_string(r.U) + ", ALPHA = " + std::to_string(r.alpha));
    const std::string output = (temporary.path() / ("itime" + std::to_string(r.U) + ".h5")).string();
    alps::solvers::ctint(parameters(r.U, r.alpha, true), output);
    if (!is_master())
      continue;
    expect_exact(scalar_result(output, "Sz2_0"), exact_sz2(r.U), "Sz2_0");
    expect_exact(scalar_result(output, "Sz0_Sz0"), exact_sz2(r.U), "Sz0_Sz0");
    expect_exact(scalar_result(output, "density_0"), -0.5, "density_0 = G(0+)");
  }
}
