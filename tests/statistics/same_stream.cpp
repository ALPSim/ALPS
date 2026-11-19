// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// One stream through every library. Quantities that are the same arithmetic
// must agree to rounding; the bins and the binned errors must equal the values
// computed by hand from the same stream. The error a library reports by default
// is not compared across libraries: they choose different bin lengths by design
// (see ar1_known_answers.cpp).
#include "streams.hpp"

#include <alps/alea.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>
#include <vector>

namespace {

using statistics_test::ar1_series;
using statistics_test::binned_error;

const std::size_t n = std::size_t(1) << 16;
const std::uint64_t seed = 20261116;

double sequential_mean(const std::vector<double> & values) {
    double sum = 0;
    for (double value : values)
        sum += value;
    return sum / values.size();
}

void expect_relative(double expected, double actual, double tolerance, const char * what, std::size_t level = 0) {
    EXPECT_NEAR(expected, actual, tolerance * std::abs(expected)) << what << ", level " << level;
}

class SameStream : public ::testing::TestWithParam<double> {
  protected:
    SameStream() : values(ar1_series(seed, GetParam(), n)), mean(sequential_mean(values)) {}
    std::vector<double> values;
    double mean;
};

}

TEST_P(SameStream, Alea)
{
    using namespace alps::alea;
    var_acc<double> simple(1);
    autocorr_acc<double> binned(1);
    batch_acc<double> batch(1);
    for (double value : values) {
        simple << value;
        binned << value;
        batch << value;
    }
    const var_result<double> simple_result = simple.finalize();
    const autocorr_result<double> binned_result = binned.finalize();
    const batch_result<double> batch_result_ = batch.finalize();
    EXPECT_EQ(simple_result.count(), n);
    EXPECT_EQ(binned_result.count(), n);
    EXPECT_EQ(batch_result_.count(), n);
    EXPECT_EQ(simple_result.mean()[0], mean);
    EXPECT_EQ(binned_result.mean()[0], mean);
    // The batches are summed in a different order.
    EXPECT_NEAR(batch_result_.mean()[0], mean, 1e-13);
    expect_relative(binned_error(values, 1), simple_result.stderror()[0], 1e-11, "naive error");

    // Every level with at least two batches is the textbook binned error.
    ASSERT_GE(binned_result.nlevel(), 17u);
    for (std::size_t level = 0; level < 16; ++level) {
        ASSERT_EQ(binned_result.batch_size(level), std::size_t(1) << level);
        expect_relative(binned_error(values, std::size_t(1) << level), binned_result.level(level).stderror()[0], 1e-10,
                        "binned error", level);
    }

    // 256 batches of 256 consecutive samples each.
    const batch_data<double> & batches = batch_result_.store();
    ASSERT_EQ(batches.num_batches(), 256u);
    EXPECT_EQ(batches.batch_size(), 256u);
    for (std::size_t i = 0; i < 256; ++i) {
        ASSERT_EQ(batches.count()(i), 256u);
        const double sum = std::accumulate(values.begin() + 256 * i, values.begin() + 256 * (i + 1), 0.);
        ASSERT_NEAR(batches.batch()(0, i), sum, 1e-11) << "batch " << i;
    }
    expect_relative(binned_error(values, 256), batch_result_.stderror()[0], 1e-10, "batch error");
}

// In the middle of a run the batches are still consecutive and of one length,
// followed by the samples that do not fill a batch yet. A second run is
// appended behind the first, and no sample is lost.
TEST(Batches, AleaKeepsTimeOrderAndAppendsRuns)
{
    using namespace alps::alea;
    const std::vector<double> values = ar1_series(seed, 0.5, 1300);
    batch_acc<double> first(1, 8), second(1, 8);
    for (std::size_t t = 0; t < 1000; ++t)
        first << values[t];
    for (std::size_t t = 1000; t < 1300; ++t)
        second << values[t];

    // 1000 samples in 8 batches: 7 of 128 and one of 104.
    const batch_result<double> run = first.result();
    EXPECT_EQ(run.store().batch_size(), 128u);
    for (std::size_t i = 0; i < 8; ++i) {
        const std::size_t begin = 128 * i, end = std::min<std::size_t>(begin + 128, 1000);
        EXPECT_EQ(run.store().count()(i), end - begin);
        EXPECT_NEAR(run.store().batch()(0, i), std::accumulate(values.begin() + begin, values.begin() + end, 0.), 1e-12);
    }

    // The second run has 4 batches of 64 and one of 44. Appended at the batch
    // length of the first run, 128, there are ten batches for eight places
    // (seven and two, and the left-over 104 and 44 samples pooled), so that
    // neighbours are added: four batches of 256 and one of 276.
    first << second.result();
    const batch_result<double> merged = first.finalize();
    EXPECT_EQ(merged.count(), 1300u);
    EXPECT_NEAR(merged.mean()[0], sequential_mean(values), 1e-14);
    EXPECT_EQ(merged.store().batch_size(), 256u);
    const std::vector<std::uint64_t> counts = {256, 256, 256, 256, 276, 0, 0, 0};
    for (std::size_t i = 0; i < 8; ++i)
        EXPECT_EQ(merged.store().count()(i), counts[i]) << "batch " << i;
    // The fourth holds the last full batch of the first run and the first 128
    // samples of the second.
    EXPECT_NEAR(merged.store().batch()(0, 3),
                std::accumulate(values.begin() + 768, values.begin() + 896, 0.)
                    + std::accumulate(values.begin() + 1000, values.begin() + 1128, 0.),
                1e-12);
}

INSTANTIATE_TEST_SUITE_P(Streams, SameStream, ::testing::Values(0.0, 0.9));
