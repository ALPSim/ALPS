// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The interaction-expansion (CT-INT) solver against exact diagonalization of
// the impurity in dmft_exact_impurity.hpp, away from half filling: the
// densities and <n_up n_down> agree with the exact ones in units of their
// error bars. The solver writes no errors for G(tau), which is not compared.
//
// The solver takes the bare G0(i w_n) = 1/(i w_n + mu0 - Delta(i w_n)) and the
// interaction U (n_up - 1/2)(n_down - 1/2), so its chemical potential in the
// exact Hamiltonian, with U n_up n_down, is mu0 + U/2. With MU = mu0, the
// impurity level EPS_f = 0 and EPSSQ_f = sum_k V^2 the solver gives G0 its
// 1/(i w_n)^2 and 1/(i w_n)^3 tails, -mu0 and mu0^2 + sum_k V^2.

#include "dmft_exact_impurity.hpp"

#include <alps/testing/temporary_directory.hpp>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <alps/solvers.hpp>

#include <complex>
#include <string>
#include <vector>

using namespace exact_impurity;

namespace {
const double mu0 = 1.;
const int n_matsubara = 500;

std::string write_g0(std::filesystem::path const & dir) {
  const std::string file = (dir / "g0.h5").string();
  std::vector<std::complex<double> > g(n_matsubara);
  for (int n = 0; n < n_matsubara; ++n) {
    const std::complex<double> iw(0., (2 * n + 1) * std::acos(-1.) / beta);
    std::complex<double> hybridization = 0.;
    for (double eps : bath) hybridization += V * V / (iw - eps);
    g[n] = 1. / (iw + mu0 - hybridization);
  }
  alps::hdf5::archive ar(file, "w");
  ar << alps::make_pvp("/G0_0", g) << alps::make_pvp("/G0_1", g);
  return file;
}
}

TEST(DmftRegression, ctint_exact_diagonalization) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  alps::params p;
  p.supply("BETA", std::to_string(beta));
  p.supply("U", std::to_string(U));
  p.supply("ALPHA", "-0.01");
  p.supply("INFILE", write_g0(temporary.path()));
  p.supply("MU", std::to_string(mu0));
  for (int f = 0; f < 2; ++f) {
    p.supply("EPS_" + std::to_string(f), "0");
    p.supply("EPSSQ_" + std::to_string(f), std::to_string(n_bath * V * V));
  }
  p.supply("N_TAU", "200");
  p.supply("N_MATSUBARA", std::to_string(n_matsubara));
  p.supply("NMATSUBARA_MEASUREMENTS", "20");
  p.supply("SWEEPS", "10000000");
  p.supply("THERMALIZATION", "1000");
  p.supply("MAX_TIME", "240");  // below the CTest TIMEOUT of 300 s
  p.supply("SEED", "13");
  const std::string output = (temporary.path() / "ctint.out.h5").string();
  alps::solvers::ctint(p, output);

  alps::hdf5::archive ar(output, "r");
  std::vector<double> densities, density_errors;
  double double_occupancy, double_occupancy_error;
  ar >> alps::make_pvp("/simulation/results/densities/mean/value", densities)
     >> alps::make_pvp("/simulation/results/densities/mean/error", density_errors)
     >> alps::make_pvp("/simulation/results/density_correlation/mean/value", double_occupancy)
     >> alps::make_pvp("/simulation/results/density_correlation/mean/error", double_occupancy_error);
  ASSERT_EQ(densities.size(), 2u);
  ASSERT_EQ(density_errors.size(), 2u);

  const exact_solution exact(mu0 + U / 2);
  mc_hypothesis::family family(3);
  for (int f = 0; f < 2; ++f)
    family.expect("density " + std::to_string(f), densities[f], density_errors[f], exact.density);
  family.expect("<n_up n_down>", double_occupancy, double_occupancy_error, exact.double_occupancy);
}
