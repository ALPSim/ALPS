/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2011 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/exchange.h>
#include <gtest/gtest.h>
#include <boost/random.hpp>
#include <algorithm>
#include <vector>

struct WeightModel {
  using weight_parameter_type = double;
  static double log_weight(double weight, double beta) { return beta * weight; }
};
TEST(ExchangeOptimization, PreservesBoundsAndHistoricalOptimizedLadder) {
  boost::mt19937 engine;
  boost::variate_generator<boost::mt19937 &, boost::uniform_real<>> random(engine,
                                                                           boost::uniform_real<>());
  alps::Parameters params;
  params["BETA_MAX"] = 5;
  params["BETA_MIN"] = 1;
  params["NUM_REPLICAS"] = 10;
  alps::parapack::exmc::inverse_temperature_set ladder(params);
  std::vector<double> original(10), weights(10);
  for (int i = 0; i < 10; ++i) {
    original[i] = ladder[i];
    weights[i] = random();
  }
  std::sort(weights.begin(), weights.end());
  for (int i = 0; i < 10; ++i)
    EXPECT_NEAR(original[i], 1 + i * 4. / 9, 1e-14);
  EXPECT_NEAR(ladder.interpolate<WeightModel>(original, weights, 2.33), .3075, 5e-5);
  EXPECT_NEAR(ladder.interpolate<WeightModel>(original, weights, 3.52), .8283, 5e-5);
  ladder.optimize_h1999<WeightModel>(weights);
  // Four-significant-digit references from the established optimization
  // fixture, independent of stream formatting and native floating precision.
  const double expected[] = {1, 1.670, 2.097, 2.422, 2.641, 2.884, 3.173, 3.780, 4.396, 5};
  const double expected_cost[] = {.03476, .03548, .03599, .03496, .03477,
                                  .03441, .03546, .03552, .03515};
  EXPECT_DOUBLE_EQ(ladder[0], 1);
  EXPECT_DOUBLE_EQ(ladder[9], 5);
  for (int i = 0; i < 10; ++i) {
    SCOPED_TRACE(i);
    EXPECT_NEAR(ladder[i], expected[i], 5e-4);
    if (i == 0)
      continue;
    EXPECT_GT(ladder[i], ladder[i - 1]);
    const auto before = ladder.interpolate<WeightModel>(original, weights, ladder[i - 1]);
    const auto after = ladder.interpolate<WeightModel>(original, weights, ladder[i]);
    const auto cost = (WeightModel::log_weight(before, ladder[i - 1]) +
                       WeightModel::log_weight(after, ladder[i])) -
                      (WeightModel::log_weight(after, ladder[i - 1]) +
                       WeightModel::log_weight(before, ladder[i]));
    EXPECT_NEAR(cost, expected_cost[i - 1], 5e-6);
  }
}
