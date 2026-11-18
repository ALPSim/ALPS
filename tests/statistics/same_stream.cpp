// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// One stream through every library. Quantities that are the same arithmetic
// must agree to rounding; the bins and the binned errors must equal the values
// computed by hand from the same stream. The error a library reports by default
// is not compared across libraries: they choose different bin lengths by design
// (see ar1_known_answers.cpp).
#include "streams.hpp"

#include <alps/alea.h>
#include <alps/alea.hpp>
#include <alps/accumulators.hpp>

#include <gtest/gtest.h>

#include <cmath>
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
    alps::SimpleRealObservable simple("simple");
    alps::RealObservable binned("binned");
    for (double value : values) {
        simple << value;
        binned << value;
    }
    EXPECT_EQ(simple.count(), n);
    EXPECT_EQ(binned.count(), n);
    EXPECT_EQ(simple.mean(), mean);
    EXPECT_EQ(binned.mean(), mean);
    expect_relative(binned_error(values, 1), simple.error(), 1e-11, "naive error");

    // ALEA's binned errors carry a factor sqrt(N / (N - 1)) on top of the
    // textbook value at every level.
    const double convention = std::sqrt(double(n) / (n - 1));
    for (std::size_t level = 0; (n >> level) >= 256; ++level)
        expect_relative(convention * binned_error(values, std::size_t(1) << level), binned.error(level), 1e-10,
                        "binned error", level);

    // The default error is the one from 128 to 255 bins: 2^16 samples make 17
    // levels, of which the top 7 are not used.
    EXPECT_EQ(binned.error(), binned.error(9));

    // 128 time-ordered bins of 512 samples.
    ASSERT_EQ(binned.bin_number(), 128u);
    ASSERT_EQ(binned.bin_size(), 512u);
    for (std::size_t b = 0; b < 128; ++b) {
        double sum = 0;
        for (std::size_t i = 0; i < 512; ++i)
            sum += values[512 * b + i];
        ASSERT_NEAR(sum, binned.bin_value(b), 1e-11) << "bin " << b;
    }
}

TEST_P(SameStream, Accumulators)
{
    namespace acc = alps::accumulators;
    acc::accumulator_set set;
    set << acc::NoBinningAccumulator<double>("simple") << acc::LogBinningAccumulator<double>("log")
        << acc::FullBinningAccumulator<double>("full");
    for (double value : values) {
        set["simple"] << value;
        set["log"] << value;
        set["full"] << value;
    }
    const acc::result_set results(set);
    for (const char * name : {"simple", "log", "full"}) {
        EXPECT_EQ(results[name].count(), n) << name;
        EXPECT_NEAR(results[name].mean<double>(), mean, 1e-13) << name;
    }
    expect_relative(binned_error(values, 1), results["simple"].error<double>(), 1e-11, "naive error");
    // The default error is the textbook value from 128 to 255 bins: 256
    // samples per bin here.
    expect_relative(binned_error(values, 512), results["log"].error<double>(), 1e-10, "binned error");
}

TEST_P(SameStream, AleaV2)
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
}

INSTANTIATE_TEST_SUITE_P(Streams, SameStream, ::testing::Values(0.0, 0.9));
