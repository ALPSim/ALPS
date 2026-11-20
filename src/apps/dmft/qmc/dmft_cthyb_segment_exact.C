// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// CT-HYB-SEGMENT against exact diagonalization of the impurity in
// dmft_exact_impurity.hpp:
//
// - G(tau) and the density at mu = U, i.e. MU = U/2 relative to half filling,
//   agree with the exact ones in units of their error bars;
// - the same model given as the DMFT loop gives it (interaction matrix,
//   absolute orbital energies and one Delta per flavor, all in HDF5) gives
//   the same Monte Carlo run;
// - with a retarded interaction K(tau), MU = -2K'(0) is half filling, as
//   particle-hole symmetry requires when mu is shifted by +K'(0).

#include "dmft_exact_impurity.hpp"

#include <alps/testing/temporary_directory.hpp>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <alps/solvers.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <string>
#include <vector>

using namespace exact_impurity;

namespace {
const int n_tau = 200;

alps::params parameters(std::string const & delta_file, double mu) {
  alps::params p;
  p.supply("BETA", std::to_string(beta));
  p.supply("U", std::to_string(U));
  p.supply("MU", std::to_string(mu));
  p.supply("FLAVORS", "2");
  p.supply("N", std::to_string(n_tau));
  p.supply("NMATSUBARA", "20");
  p.supply("MAX_TIME", "600");
  p.supply("SEED", "7");
  p.supply("cthyb.DELTA", delta_file);
  p.supply("cthyb.N_MEAS", "20");
  p.supply("cthyb.SWEEPS", "4000000");
  p.supply("cthyb.THERMALIZATION", "20000");
  p.supply("cthyb.MEASURE_freq", "0");
  return p;
}

struct result {
  std::vector<double> green, green_error;
  double density[2], density_error[2];
};

// /G_tau holds the means of the measured g_0 with the end points replaced by
// the densities; its errors are those of g_0.
result read_result(std::string const & file) {
  alps::hdf5::archive ar(file, "r");
  result r;
  ar >> alps::make_pvp("/G_tau/0/mean/value", r.green)
     >> alps::make_pvp("/simulation/results/g_0/mean/error", r.green_error);
  for (int f = 0; f < 2; ++f) {
    const std::string name = "/simulation/results/density_" + std::to_string(f) + "/mean/";
    ar >> alps::make_pvp(name + "value", r.density[f]) >> alps::make_pvp(name + "error", r.density_error[f]);
  }
  return r;
}
}

TEST(DmftRegression, cthyb_segment_exact_diagonalization) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  const std::string output = (temporary.path() / "plain.out.h5").string();
  // MU is relative to half filling: the solver adds U/2.
  alps::solvers::cthyb_segment(parameters(write_delta_text(temporary.path(), n_tau), 1.), output);
  if (!is_master()) return;

  const exact_solution exact(1. + U / 2);
  const result r = read_result(output);
  ASSERT_EQ(r.green.size(), static_cast<std::size_t>(n_tau + 1));
  ASSERT_EQ(r.green_error.size(), r.green.size());
  mc_hypothesis::family family(2 + 9);
  for (int f = 0; f < 2; ++f)
    family.expect("density_" + std::to_string(f), r.density[f], r.density_error[f], exact.density);
  for (int t = 20; t < n_tau; t += 20)
    family.expect("G(tau = " + std::to_string(beta * t / n_tau) + ")", r.green[t], r.green_error[t],
                  exact.green_bin(t, n_tau));
  EXPECT_NEAR(r.green[0], -(1 - r.density[0]), 1e-12);
  EXPECT_NEAR(r.green[n_tau], -r.density[0], 1e-12);
}

TEST(DmftRegression, cthyb_segment_dmft_inputs) {
#ifdef ALPS_HAVE_MPI
  // Under MPI the run stops at scheduled checks, not after exactly SWEEPS.
  GTEST_SKIP() << "compares two runs update by update, which needs a serial build";
#endif
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  const std::string plain = (temporary.path() / "plain.out.h5").string();
  alps::solvers::cthyb_segment(parameters(write_delta_text(temporary.path(), n_tau), 1.), plain);

  // The interaction matrix, absolute orbital energies and Delta per flavor
  // in one HDF5 file, as the DMFT loop and the tutorials provide them.
  // Every rank has its own temporary directory and writes its own inputs.
  const std::string input = (temporary.path() / "input.h5").string();
  {
    alps::hdf5::archive ar(input, "w");
    std::vector<double> d(n_tau + 1);
    for (int i = 0; i <= n_tau; ++i) d[i] = delta(beta * i / n_tau);
    ar << alps::make_pvp("/Delta_0", d) << alps::make_pvp("/Delta_1", d);
    ar << alps::make_pvp("/Umatrix", std::vector<double>{0., U, U, 0.});
    ar << alps::make_pvp("/MUvector", std::vector<double>{1. + U / 2, 1. + U / 2});
  }
  alps::params p = parameters(input, 0.);
  p.supply("cthyb.DELTA_IN_HDF5", "1");
  p.supply("U_MATRIX", input);
  p.supply("UMATRIX_IN_HDF5", "1");
  p.supply("MU_VECTOR", input);
  p.supply("MU_IN_HDF5", "1");
  const std::string from_hdf5 = (temporary.path() / "hdf5.out.h5").string();
  alps::solvers::cthyb_segment(p, from_hdf5);
  if (!is_master()) return;

  const result a = read_result(plain), b = read_result(from_hdf5);
  for (int f = 0; f < 2; ++f) EXPECT_EQ(a.density[f], b.density[f]);
  EXPECT_EQ(a.green, b.green);
}

TEST(DmftRegression, cthyb_segment_retarded_interaction_half_filling) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  // K(tau) of a Holstein phonon of frequency w0 and coupling lambda; K'(0) = lambda^2/w0
  const double lambda = 1., w0 = 5.;
  const std::string kfile = (temporary.path() / "K.h5").string();
  {
    std::vector<double> K(n_tau + 1), Kp(n_tau + 1);
    for (int i = 0; i <= n_tau; ++i) {
      const double tau = beta * i / n_tau, s = std::sinh(w0 * beta / 2);
      K[i] = -lambda * lambda / (w0 * w0) * (std::cosh(w0 * (beta / 2 - tau)) / s - std::cosh(w0 * beta / 2) / s);
      Kp[i] = lambda * lambda / w0 * std::sinh(w0 * (beta / 2 - tau)) / s;
    }
    alps::hdf5::archive ar(kfile, "w");
    ar << alps::make_pvp("/Ret_int_K", K) << alps::make_pvp("/Ret_int_Kp", Kp);
  }
  // With mu -> mu-K'(0) instead of mu+K'(0), the density is 0.44.
  alps::params p = parameters(write_delta_text(temporary.path(), n_tau), -2 * lambda * lambda / w0);
  p.supply("cthyb.RET_INT_K", kfile);
  p.supply("cthyb.K_IN_HDF5", "1");
  const std::string output = (temporary.path() / "retarded.out.h5").string();
  alps::solvers::cthyb_segment(p, output);
  if (!is_master()) return;

  const result r = read_result(output);
  // The run polarizes the spin slowly, more slowly than the binning analysis
  // sees, so the two densities deviate in opposite directions by more than
  // their errors. Their average cancels that; the mean of the two errors
  // bounds its error from above whatever their correlation.
  mc_hypothesis::family family(1);
  family.expect("(density_0 + density_1) / 2", (r.density[0] + r.density[1]) / 2,
                (r.density_error[0] + r.density_error[1]) / 2, 0.5);
}
