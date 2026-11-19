// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// The definitions, on a sample whose sums are small integers: count, mean,
// variance and the naive standard error of the mean. Every expected value is
// exact, so a change of normalisation in any library fails here first.
#include <alps/alea.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

namespace {

// x = 1..8: mean 4.5, sum of squared deviations 42, unbiased variance 6.
const int n = 8;
const double mean = 4.5;
const double variance = 6.0;
const double naive_error = std::sqrt(6.0 / 8.0);

}

TEST(ExactMoments, Alea)
{
    using namespace alps::alea;
    var_acc<double> simple(1);
    autocorr_acc<double> binned(1);
    batch_acc<double> batch(1);
    var_acc<double> vector(2);
    for (int t = 1; t <= n; ++t) {
        simple << double(t);
        binned << double(t);
        batch << double(t);
        vector << std::vector<double>({double(t), 2. * t});
    }
    const var_result<double> simple_result = simple.finalize();
    EXPECT_EQ(simple_result.count(), 8u);
    EXPECT_EQ(simple_result.mean()[0], mean);
    EXPECT_DOUBLE_EQ(simple_result.var()[0], variance);
    EXPECT_DOUBLE_EQ(simple_result.stderror()[0], naive_error);

    const autocorr_result<double> binned_result = binned.finalize();
    EXPECT_EQ(binned_result.count(), 8u);
    EXPECT_EQ(binned_result.mean()[0], mean);
    EXPECT_DOUBLE_EQ(binned_result.level(0).stderror()[0], naive_error);
    // Level 1 holds the four pair means 1.5, 3.5, 5.5, 7.5: sum of squared
    // deviations 20, standard error sqrt(20 / (4 * 3)).
    EXPECT_DOUBLE_EQ(binned_result.level(1).stderror()[0], std::sqrt(20.0 / 12.0));

    const batch_result<double> batch_result_ = batch.finalize();
    EXPECT_EQ(batch_result_.count(), 8u);
    EXPECT_EQ(batch_result_.mean()[0], mean);

    const var_result<double> vector_result = vector.finalize();
    EXPECT_EQ(vector_result.mean()[1], 2 * mean);
    EXPECT_DOUBLE_EQ(vector_result.stderror()[1], 2 * naive_error);
}
