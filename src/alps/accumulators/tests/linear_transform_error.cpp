/*
 * Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
 */

// A full-binning result shifted, scaled or negated by a constant keeps the
// binning-analysis error of the result, scaled by the constant's magnitude.
// The data are correlated, so a jackknife over the bins would give another
// error.

#include "alps/accumulators.hpp"
#include "gtest/gtest.h"

#include <random>
#include <utility>
#include <vector>

namespace acc = alps::accumulators;

namespace {

acc::result_set correlated_results() {
    acc::accumulator_set measurements;
    measurements << acc::FullBinningAccumulator<double>("x")
                 << acc::FullBinningAccumulator<std::vector<double> >("v");
    std::mt19937 engine(1);
    std::normal_distribution<double> xi;
    double x = 0;
    for (int i = 0; i < 2000; ++i) {
        x = 0.8 * x + xi(engine);
        measurements["x"] << x;
        measurements["v"] << std::vector<double>({x, 2 * x});
    }
    return acc::result_set(measurements);
}

} // namespace

TEST(LinearTransformError, ScalarKeepsTheError) {
    acc::result_set results = correlated_results();
    acc::result_wrapper const& r = results["x"];
    double const mean = r.mean<double>(), error = r.error<double>();

    EXPECT_DOUBLE_EQ(error, (-r).error<double>());
    EXPECT_DOUBLE_EQ(error, (r + 1.0).error<double>());
    EXPECT_DOUBLE_EQ(error, (r - 1.0).error<double>());
    EXPECT_DOUBLE_EQ(2 * error, (r * 2.0).error<double>());
    EXPECT_DOUBLE_EQ(2 * error, (r * -2.0).error<double>());
    EXPECT_DOUBLE_EQ(error / 2, (r / 2.0).error<double>());
    EXPECT_DOUBLE_EQ(error, (r * 2.0 / 2.0).error<double>());
    EXPECT_DOUBLE_EQ(error, (1.0 - r).error<double>());
    EXPECT_DOUBLE_EQ(2 * error, (2.0 * r).error<double>());

    EXPECT_DOUBLE_EQ(-mean, (-r).mean<double>());
    EXPECT_DOUBLE_EQ(mean + 1, (r + 1.0).mean<double>());
    EXPECT_DOUBLE_EQ(2 * mean, (r * 2.0).mean<double>());
}

TEST(LinearTransformError, VectorKeepsTheError) {
    acc::result_set results = correlated_results();
    acc::result_wrapper const& r = results["v"];
    std::vector<double> const error = r.error<std::vector<double> >();

    std::vector<double> const negated = (-r).error<std::vector<double> >();
    std::vector<double> const scaled = (r * 3.0).error<std::vector<double> >();
    ASSERT_EQ(2u, negated.size());
    ASSERT_EQ(2u, scaled.size());
    for (std::size_t i = 0; i < 2; ++i) {
        EXPECT_DOUBLE_EQ(error[i], negated[i]);
        EXPECT_DOUBLE_EQ(3 * error[i], scaled[i]);
    }
}

// The count and mean keep the samples of the last, incomplete bin.
TEST(LinearTransformError, IncompleteBinIsKept) {
    acc::accumulator_set measurements;
    measurements << acc::FullBinningAccumulator<double>("x");
    for (int i = 0; i < 129; ++i)
        measurements["x"] << double(i);
    acc::result_set results(measurements);
    acc::result_wrapper const& r = results["x"];
    ASSERT_EQ(129u, r.count());
    ASSERT_DOUBLE_EQ(64.0, r.mean<double>());
    std::vector<std::pair<acc::result_wrapper, double> > const transformed = {
        {-r, -64.0}, {r + 1.0, 65.0}, {r * 2.0, 128.0}, {1.0 - r, -63.0}};
    for (auto const& t : transformed) {
        EXPECT_EQ(129u, t.first.count());
        EXPECT_DOUBLE_EQ(t.second, t.first.mean<double>());
    }
}

// A nonlinear transform still uses the jackknife, also after a linear one.
TEST(LinearTransformError, NonlinearAfterLinearUsesTheJackknife) {
    acc::result_set results = correlated_results();
    acc::result_wrapper const& r = results["x"];
    acc::result_wrapper const squared = r * r;
    acc::result_wrapper const scaled_squared = (r * 2.0) * (r * 2.0);
    EXPECT_NEAR(4 * squared.mean<double>(), scaled_squared.mean<double>(), 1e-12);
    EXPECT_NEAR(4 * squared.error<double>(), scaled_squared.error<double>(), 1e-12);
}
