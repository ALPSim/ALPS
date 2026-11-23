// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The hybridization (CT-HYB) solver against exact diagonalization of the
// impurity in dmft_exact_impurity.hpp, away from half filling: the densities
// and G(tau) at nine imaginary times agree with the exact ones in units of
// their error bars. The sweeps make an off-by-one in the segment insertion's
// proposal ratio about ten error bars in the densities; it runs about 8 s in
// a Release build, the solver being slower per sweep than CT-HYB-SEGMENT.
// Unoptimized builds are about 30 times slower and run a third of the sweeps;
// the comparison stays valid with larger error bars.

#include "dmft_exact_impurity.hpp"

#include <alps/testing/temporary_directory.hpp>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <alps/solvers.hpp>

#include <string>
#include <vector>

using namespace exact_impurity;

namespace {
const int n_tau = 200;
const double mu = 1. + U / 2;
#ifdef NDEBUG
const char * const sweeps = "600000";
#else
const char * const sweeps = "200000";
#endif
}

TEST(DmftRegression, cthyb_exact_diagonalization) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  alps::params p;
  p.supply("BETA", std::to_string(beta));
  p.supply("U", std::to_string(U));
  p.supply("MU", std::to_string(mu));  // the hybridization solver's MU is absolute
  p.supply("N_ORBITALS", "2");
  p.supply("N_TAU", std::to_string(n_tau));
  p.supply("DELTA", write_delta_text(temporary.path(), n_tau));
  p.supply("N_MEAS", "20");
  p.supply("SWEEPS", sweeps);
  p.supply("THERMALIZATION", "20000");
  p.supply("MAX_TIME", "240");  // below the CTest TIMEOUT of 300 s
  p.supply("SEED", "11");
  const std::string output = (temporary.path() / "cthyb.out.h5").string();
  alps::solvers::cthyb(p, output);
  if (!is_master()) return;

  alps::hdf5::archive ar(output, "r");
  std::vector<double> green, green_error;
  ar >> alps::make_pvp("/G_tau/0/mean/value", green) >> alps::make_pvp("/G_tau/0/mean/error", green_error);
  ASSERT_EQ(green.size(), static_cast<std::size_t>(n_tau + 1));
  ASSERT_EQ(green_error.size(), green.size());

  const exact_solution exact(mu);
  mc_hypothesis::family family(2 + 9);
  for (int f = 0; f < 2; ++f) {
    const std::string name = "/simulation/results/density_" + std::to_string(f) + "/mean/";
    double density, error;
    ar >> alps::make_pvp(name + "value", density) >> alps::make_pvp(name + "error", error);
    family.expect("density_" + std::to_string(f), density, error, exact.density);
  }
  for (int t = 20; t < n_tau; t += 20)
    family.expect("G(tau = " + std::to_string(beta * t / n_tau) + ")", green[t], green_error[t],
                  exact.green_bin(t, n_tau));
}
