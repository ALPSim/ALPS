// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The interaction-expansion solver measures Sz_i, Sz2_i and Sz0_Szi as the
// value times the sign. ALPS 3.0 measured them as signed observables, so
// /simulation/results/Sz_0/mean/value held <Sz s>/<s>. The solver divides
// the results by "Sign" with divide_sz_results_by_sign before saving them;
// this test pins that with a mean sign of 1/2.

#include "interaction_expansion2/sz_results.hpp"

#include <gtest/gtest.h>
#include <alps/hdf5/archive.hpp>

#include <cstdio>
#include <string>

TEST(DmftRegression, ctint_sz_results_divided_by_sign) {
  using alps::accumulators::FullBinningAccumulator;
  alps::accumulators::accumulator_set measurements;
  for (char const * name : {"Sign", "Sz_0", "Sz2_0", "Sz0_Sz0", "Sz_1", "density_0"})
    measurements << FullBinningAccumulator<double>(name);

  // The samples repeat with period 4, which divides every bin, so the
  // jackknife bins are equal and the divided means are exact ratios.
  double const sign[4] = {1., 1., 1., -1.};
  double const sz[4] = {0.5, -0.25, 1., 0.75};
  int const n = 4096;
  for (int k = 0; k < n; ++k) {
    double const s = sign[k % 4], x = sz[k % 4];
    measurements["Sign"] << s;
    measurements["Sz_0"] << x * s;
    measurements["Sz2_0"] << x * x * s;
    measurements["Sz0_Sz0"] << x * x * s;
    measurements["density_0"] << x * s;
  }

  alps::accumulators::result_set results(measurements);
  ASSERT_EQ(results["Sz_1"].count(), 0u);
  divide_sz_results_by_sign(results, 2);

  // <s> = 1/2, <Sz s> = 1/8 and <Sz^2 s> = 3/16 over one period.
  EXPECT_DOUBLE_EQ(results["Sign"].mean<double>(), 0.5);
  EXPECT_DOUBLE_EQ(results["Sz_0"].mean<double>(), 0.25);
  EXPECT_DOUBLE_EQ(results["Sz2_0"].mean<double>(), 0.375);
  EXPECT_DOUBLE_EQ(results["Sz0_Sz0"].mean<double>(), 0.375);
  EXPECT_EQ(results["Sz_0"].count(), static_cast<unsigned>(n));
  // Unmeasured Sz results and other observables are left alone.
  EXPECT_EQ(results["Sz_1"].count(), 0u);
  EXPECT_DOUBLE_EQ(results["density_0"].mean<double>(), 0.125);

  std::string const file = "ctint_sz_results.h5";
  {
    alps::hdf5::archive ar(file, "w");
    ar["/simulation/results"] << results;
  }
  double saved = 0.;
  bool sign_changed = true;
  {
    alps::hdf5::archive ar(file, "r");
    ar["/simulation/results/Sz_0/mean/value"] >> saved;
    sign_changed = ar.is_group("/simulation/results/Sign/jacknife");
  }
  std::remove(file.c_str());
  EXPECT_DOUBLE_EQ(saved, 0.25);
  EXPECT_FALSE(sign_changed);
}
