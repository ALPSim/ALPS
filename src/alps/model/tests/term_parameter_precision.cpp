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

// Parameter values substituted into SITETERM and BONDTERM expressions keep
// full double precision. They used to be printed back into the term with 6
// significant digits, so mu=0.123456789123 acted as 0.123457.

#include <gtest/gtest.h>
#include <alps/lattice.h>
#include <alps/model.h>
#include <alps/model/hamiltonian_matrix.hpp>
#include <alps/numeric/matrix.hpp>

#include <cmath>

TEST(TermParameterPrecision, SiteAndBondTermsKeepFullPrecision) {
  const double mu = 0.123456789123, t = 0.987654321987;
  alps::Parameters p;
  p["LATTICE"] = "open chain lattice";
  p["L"] = 2;
  p["MODEL"] = "fermion Hubbard";
  p["CONSERVED_QUANTUMNUMBERS"] = "Nup,Ndown";
  p["Nup_total"] = 1;
  p["Ndown_total"] = 0;
  p["mu"] = mu;
  p["t"] = t;
  // One up electron on two sites: H = [[-mu, -/+t], [-/+t, -mu]].
  alps::hamiltonian_matrix<alps::numeric::matrix<double> > h(p);
  const alps::numeric::matrix<double>& m = h.matrix();
  ASSERT_EQ(h.dimension(), 2u);
  EXPECT_NEAR(m(0, 0), -mu, 1e-15);
  EXPECT_NEAR(m(1, 1), -mu, 1e-15);
  EXPECT_NEAR(std::abs(m(0, 1)), t, 1e-15);
  EXPECT_NEAR(std::abs(m(1, 0)), t, 1e-15);
}
