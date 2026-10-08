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

TEST(LinearRegression, HistoricalSubrangesAgreeWithCenteredReference) {
  boost::mt19937 engine;
  boost::variate_generator<boost::mt19937 &, boost::uniform_real<>> random(engine,
                                                                           boost::uniform_real<>());
  std::vector<double> x, y;
  for (int i = 0; i < 10; ++i) {
    x.push_back(random());
    y.push_back(random());
  }
  std::sort(x.begin(), x.end());
  std::sort(y.begin(), y.end());
  for (const auto selection : {std::pair<int, int>{10, 0}, {2, 2}, {5, 3}}) {
    const int count = selection.first, start = selection.second;
    SCOPED_TRACE(start);
    SCOPED_TRACE(count);
    long double mean_x = 0, mean_y = 0;
    for (int i = start; i < start + count; ++i) {
      mean_x += x[i];
      mean_y += y[i];
    }
    mean_x /= count;
    mean_y /= count;
    long double covariance = 0, variance = 0;
    for (int i = start; i < start + count; ++i) {
      covariance += (x[i] - mean_x) * (y[i] - mean_y);
      variance += (x[i] - mean_x) * (x[i] - mean_x);
    }
    const auto result =
        alps::parapack::exmc::inverse_temperature_set::linear_regression(count, start, x, y);
    EXPECT_NEAR(result.first, static_cast<double>(mean_y - covariance / variance * mean_x), 1e-12);
    EXPECT_NEAR(result.second, static_cast<double>(covariance / variance), 1e-12);
  }
}
TEST(LinearRegression, RecoversAnExactLine) {
  const std::vector<double> x{0, 1, 2, 3, 4}, y{2, 5, 8, 11, 14};
  const auto result = alps::parapack::exmc::inverse_temperature_set::linear_regression(5, 0, x, y);
  EXPECT_DOUBLE_EQ(result.first, 2);
  EXPECT_DOUBLE_EQ(result.second, 3);
}
