/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// The global K of the `spin` model sets the biquadratic coupling of the
// nearest-neighbour bonds, as J does for the exchange. It used to be dropped,
// so only K0 had an effect.

#include <gtest/gtest.h>
#include <alps/lattice.h>
#include <alps/model.h>
#include <alps/model/hamiltonian_matrix.hpp>
#include <alps/numeric/matrix.hpp>

#include <cmath>
#include <cstddef>
#include <string>

namespace {
alps::numeric::matrix<double> spin1_dimer(std::string const& name) {
  alps::Parameters p;
  p["LATTICE"] = "open chain lattice";
  p["L"] = 2;
  p["MODEL"] = "spin";
  p["local_S"] = 1;
  p["J"] = 0;
  p[name] = 1;
  return alps::hamiltonian_matrix<alps::numeric::matrix<double> >(p).matrix();
}
}

TEST(SpinModel, GlobalKSetsTheBiquadraticCoupling) {
  alps::numeric::matrix<double> const k = spin1_dimer("K");
  alps::numeric::matrix<double> const k0 = spin1_dimer("K0");
  ASSERT_EQ(k.num_rows(), 9u);
  ASSERT_EQ(k0.num_rows(), 9u);
  double norm = 0;
  for (std::size_t i = 0; i < 9; ++i)
    for (std::size_t j = 0; j < 9; ++j) {
      EXPECT_NEAR(k(i, j), k0(i, j), 1e-14);
      norm += std::abs(k(i, j));
    }
  EXPECT_GT(norm, 1);
}
