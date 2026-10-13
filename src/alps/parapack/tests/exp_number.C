/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/parapack/exp_number.h>
#include <gtest/gtest.h>
#include <cmath>

TEST(ExponentialNumber, MixedArithmeticAgreesWithRealArithmetic) {
  const alps::exp_double x = 3, y = 2, p = -2.5;
  // Ordinary magnitudes use explicit absolute roundoff bounds; huge values
  // are tested in logarithmic representation in the separate case below.
  EXPECT_NEAR(static_cast<double>(x), 3, 1e-12);
  EXPECT_NEAR(static_cast<double>(y), 2, 1e-12);
  EXPECT_NEAR(static_cast<double>(x + y), 5, 1e-12);
  EXPECT_NEAR(static_cast<double>(x * y), 6, 1e-12);
  EXPECT_NEAR(static_cast<double>(x / y), 1.5, 1e-12);
  EXPECT_NEAR(static_cast<double>(x + 1.2), 4.2, 1e-12);
  EXPECT_NEAR(static_cast<double>(x - 1.2), 1.8, 1e-12);
  EXPECT_NEAR(static_cast<double>(x - 3.5), -.5, 1e-12);
  EXPECT_NEAR(static_cast<double>(x * 1.2), 3.6, 1e-12);
  EXPECT_NEAR(static_cast<double>(x * 2), 6, 1e-12);
  EXPECT_NEAR(static_cast<double>(x / 1.2), 2.5, 1e-12);
  EXPECT_NEAR(static_cast<double>(x / 3), 1, 1e-12);
  EXPECT_NEAR(static_cast<double>(3.5 + x), 6.5, 1e-12);
  EXPECT_NEAR(static_cast<double>(3 + x), 6, 1e-12);
  EXPECT_NEAR(static_cast<double>(3.5 - x), .5, 1e-12);
  EXPECT_NEAR(static_cast<double>(5 - x), 2, 1e-12);
  EXPECT_NEAR(static_cast<double>(4 - x), 1, 1e-12);
  EXPECT_NEAR(static_cast<double>(1.2 - x), -1.8, 1e-12);
  EXPECT_NEAR(static_cast<double>(3.5 * x), 10.5, 1e-12);
  EXPECT_NEAR(static_cast<double>(3 * x), 9, 1e-12);
  EXPECT_NEAR(static_cast<double>(3.5 / x), 3.5 / 3, 1e-12);
  EXPECT_NEAR(static_cast<double>(3 / x), 1, 1e-12);
  EXPECT_NEAR(static_cast<double>(p), -2.5, 1e-12);
  EXPECT_NEAR(static_cast<double>(-p), 2.5, 1e-12);
  EXPECT_EQ(static_cast<double>(p - p), 0);
  EXPECT_NEAR(static_cast<double>(x + p), .5, 1e-12);
  EXPECT_NEAR(static_cast<double>(x - p), 5.5, 1e-12);
  EXPECT_NEAR(static_cast<double>(x * p), -7.5, 1e-12);
  EXPECT_NEAR(static_cast<double>(x / p), -1.2, 1e-12);
  EXPECT_NEAR(static_cast<double>(y + p), -.5, 1e-12);
  EXPECT_NEAR(static_cast<double>(y - p), 4.5, 1e-12);
  EXPECT_NEAR(static_cast<double>(y * p), -5, 1e-12);
  EXPECT_NEAR(static_cast<double>(y / p), -.8, 1e-12);
}

TEST(ExponentialNumber, HugeMagnitudesRemainUsableInLogSpace) {
  const auto v = alps::exp_value(10000), w = alps::exp_value(9999);
  EXPECT_TRUE(std::isinf(static_cast<double>(v)));
  EXPECT_TRUE(std::isinf(static_cast<double>(w)));
  EXPECT_DOUBLE_EQ(log(pow(v, 3)), 30000);
  EXPECT_NEAR(log(v + w), 10000 + std::log1p(std::exp(-1.)), 1e-11);
  EXPECT_NEAR(log(v - w), 10000 + std::log1p(-std::exp(-1.)), 1e-11);
  EXPECT_DOUBLE_EQ(log(v * w), 19999);
  EXPECT_NEAR(static_cast<double>(v / w), std::exp(1.), 1e-12);
  EXPECT_GT(v, w);
  EXPECT_GE(v, w);
  // The legacy API supplies operator== but no operator!= in C++17.
  // Test its equality contract directly, avoiding conversion of both to inf.
  EXPECT_FALSE(v == w);
  EXPECT_FALSE(v <= w);
  EXPECT_FALSE(v < w);
  EXPECT_FALSE(v > v);
  EXPECT_GE(v, v);
  EXPECT_EQ(v, v);
  EXPECT_LE(v, v);
  EXPECT_FALSE(v < v);
}

TEST(ExponentialNumber, SignedAndMixedComparisons) {
  const alps::exp_double x = 3, y = 2, p = -2.5;
  for (const auto value : {x, y}) {
    EXPECT_GT(value, p);
    EXPECT_GE(value, p);
    EXPECT_FALSE(value == p);
    EXPECT_FALSE(value <= p);
    EXPECT_FALSE(value < p);
  }
  EXPECT_EQ(p, p);
  EXPECT_FALSE(p == -p);
  EXPECT_TRUE(x > 2);
  EXPECT_TRUE(x >= 2);
  EXPECT_FALSE(x == 2);
  EXPECT_FALSE(x <= 2);
  EXPECT_FALSE(x < 2);
}
