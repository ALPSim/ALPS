/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/measurement.h>
#include <gtest/gtest.h>

TEST(MeasurementMerge, PreservesHistoricalMeansAndBinningErrors) {
  alps::RealObservable single("single");
  for (int i = 0; i < 100; ++i)
    single << static_cast<double>(i);
  EXPECT_DOUBLE_EQ(single.mean(), 49.5);
  // Legacy output rounded errors to three significant digits. Preserve that
  // uncertainty contract without coupling the test to ostream formatting.
  EXPECT_NEAR(single.error(), 2.92, .005);
  alps::ObservableSet merged, random_clone;
  for (int i = 0; i < 100; ++i) {
    alps::ObservableSet bin;
    bin << alps::RealObservable("obs", 10000);
    bin.reset(true);
    for (int j = 0; j < 100; ++j)
      bin["obs"] << static_cast<double>(i);
    merged << bin;
    alps::ObservableSet clone;
    clone << alps::RealObservable("obs");
    clone.reset(true);
    for (int j = 0; j < 100; ++j)
      clone["obs"] << static_cast<double>(i);
    alps::merge_random_clone(random_clone, clone);
  }
  const auto &merged_result = dynamic_cast<const alps::RealObsevaluator &>(merged["obs"]);
  // The evaluator combines 10,000 samples through jackknife arithmetic;
  // allow accumulated roundoff while tightening the old 0.05 output bound.
  EXPECT_NEAR(merged_result.mean(), 49.5, 1e-8);
  EXPECT_NEAR(merged_result.error(), .289, .0005);
  const auto &clone_result = dynamic_cast<const alps::RealObservable &>(random_clone["obs"]);
  EXPECT_DOUBLE_EQ(clone_result.mean(), single.mean());
  EXPECT_DOUBLE_EQ(clone_result.error(), single.error());
}
