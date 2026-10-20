/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2012 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/random/uniform_on_sphere_n.h>
#include <alps/utilities/vectorio.hpp>
#include <boost/random.hpp>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include <numeric>

TEST(UniformOnSphere, SeededSamplesMatchReferenceAndHaveUnitNorm) {
  boost::mt19937 eng(5489);
  alps::uniform_on_sphere_n<3, double, std::vector<double>> dist;
  // Preserve the 100 historical seed-5489 samples, printed to three
  // significant digits in uniform_on_sphere_n.op.
  std::ifstream reference(std::string(ALPS_TEST_SOURCE_DIR) + "/uniform_on_sphere_n.op");
  ASSERT_TRUE(reference.good());
  for (int i = 0; i < 100; ++i) {
    const auto r = dist(eng);
    ASSERT_EQ(r.size(), 3u);
    char delimiter;
    ASSERT_TRUE(static_cast<bool>(reference >> delimiter));
    ASSERT_EQ(delimiter, '(');
    for (int axis = 0; axis < 3; ++axis) {
      double expected;
      ASSERT_TRUE(static_cast<bool>(reference >> expected >> delimiter));
      ASSERT_EQ(delimiter, axis == 2 ? ')' : ',');
      ASSERT_NE(expected, 0.0); // All historical coordinates are nonzero.
      const double tolerance = 0.500001 * std::pow(10., std::floor(std::log10(std::abs(expected))) - 2);
      EXPECT_NEAR(r[axis], expected, tolerance) << "sample=" << i << ", axis=" << axis;
    }
    // Squared norm accumulates three rounded products; 1e-12 allows roundoff.
    EXPECT_NEAR(std::inner_product(r.begin(), r.end(), r.begin(), 0.), 1., 1e-12)
      << "seed=5489, sample=" << i;
  }
  reference >> std::ws;
  EXPECT_TRUE(reference.eof());
}
