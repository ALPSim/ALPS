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

#include <alps/parapack/wanglandau.h>
#include <alps/alea.h>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>
#include <boost/random.hpp>
#include <cmath>
#include <vector>

TEST(WangLandauWeights, CountsVisitsAndMergesOverlappingSavedRanges) {
  boost::mt19937 engine;
  boost::variate_generator<boost::mt19937 &, boost::uniform_real<>> random(engine,
                                                                           boost::uniform_real<>());
  alps::testing::TemporaryDirectory directory;
  alps::Parameters params;
  params["DIR_NAME"] = directory.path().string();
  const std::vector<double> expected0{3.46574,  1.38629, 2.77259,  4.85203, 2.77259, 2.07944,
                                      1.38629,  3.46574, 0.693147, 2.77259, 3.46574, 3.46574,
                                      0.693147, 1.38629, 0.0,      4.15888, 2.77259, 2.77259,
                                      2.77259,  7.62462, 2.07944,  2.77259, 3.46574, 6.23832};
  const std::vector<double> expected1{27.0327, 31.1916, 33.2711, 31.8848, 35.3505, 38.8162,
                                      34.6574, 32.5779, 34.6574, 36.7368, 36.7368, 32.5779,
                                      33.2711, 34.6574, 32.5779, 31.8848, 31.1916, 31.1916,
                                      26.3396, 31.8848, 34.6574};
  const std::vector<double> expected2{30.6718, 28.5923, 29.9786, 32.0581, 29.9786, 29.2855, 28.5923,
                                      30.6718, 27.8992, 29.9786, 30.6718, 30.6718, 27.8992, 28.5923,
                                      0.0,     31.3649, 29.9786, 29.9786, 29.9786, 34.8306, 28.1591,
                                      30.5851, 31.9714, 32.6646, 35.3505, 38.8162, 34.6574, 32.5779,
                                      34.6574, 36.7368, 36.7368, 32.5779, 33.2711, 34.6574, 32.5779,
                                      31.8848, 31.1916, 31.1916, 26.3396, 31.8848, 34.6574};
  const alps::integer_range<int> range0(-10, 13), range1(10, 30);
  alps::wanglandau_weight weight0(range0), weight1(range1);
  alps::ObservableSet observables0, observables1;
  weight0.init_observables(observables0);
  observables0.reset(true);
  weight1.init_observables(observables1);
  observables1.reset(true);
  std::vector<int> visits0(range0.size()), visits1(range1.size());
  for (int i = 0; i < 100; ++i) {
    const int index = static_cast<int>(range0.size() * random());
    ++visits0[index];
    weight0.visit(observables0, index + range0.min(), 2);
  }
  for (int i = 0; i < 1000; ++i) {
    const int index = static_cast<int>(range1.size() * random());
    ++visits1[index];
    weight1.visit(observables1, index + range1.min(), 2);
  }
  for (int bin = range0.min(); bin <= range0.max(); ++bin) {
    SCOPED_TRACE(bin);
    EXPECT_NEAR(log(weight0[bin]), visits0[bin - range0.min()] * std::log(2.), 1e-12);
    EXPECT_NEAR(log(weight0[bin]), expected0[bin - range0.min()], 5e-5);
  }
  for (int bin = range1.min(); bin <= range1.max(); ++bin) {
    SCOPED_TRACE(bin);
    EXPECT_NEAR(log(weight1[bin]), visits1[bin - range1.min()] * std::log(2.), 1e-12);
    EXPECT_NEAR(log(weight1[bin]), expected1[bin - range1.min()], 5e-5);
  }
  weight0.write_observables(observables0);
  weight1.write_observables(observables1);
  params["WEIGHT_DUMP_FILE"] = "weight0.xdr";
  alps::wanglandau_weight::save_weight(observables0, params);
  params["WEIGHT_DUMP_FILE"] = "weight1.xdr";
  alps::wanglandau_weight::save_weight(observables1, params);
  alps::wanglandau_weight merged(unify(range0, range1));
  params["WEIGHT_DUMP_FILE"] = "weight0.xdr weight1.xdr";
  merged.load_weight(params);
  for (int bin = range0.min(); bin <= range1.max(); ++bin) {
    SCOPED_TRACE(bin);
    // Six-significant-digit historical merged values, including the unvisited
    // bin and the four-bin overlap. Normalization changes are observable here.
    EXPECT_NEAR(log(merged[bin]), expected2[bin - range0.min()], 5e-5);
  }
}
