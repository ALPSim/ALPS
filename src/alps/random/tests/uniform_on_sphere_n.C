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
#include <alps/utility/vectorio.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <vector>

#include <gtest/gtest.h>
#include <numeric>

TEST(UniformOnSphere, SamplesHaveUnitNorm) {
  boost::mt19937 eng(5489);
  alps::uniform_on_sphere_n<3, double, std::vector<double>> dist;
  for (int i = 0; i < 100; ++i) {
    const auto r = dist(eng);
    ASSERT_EQ(r.size(), 3u);
    // Squared norm accumulates three rounded products; 1e-12 allows roundoff.
    EXPECT_NEAR(std::inner_product(r.begin(), r.end(), r.begin(), 0.), 1., 1e-12)
      << "seed=5489, sample=" << i;
  }
}
