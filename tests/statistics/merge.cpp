// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Merging runs. Splitting one stream over two accumulators and merging them
// must give the count and mean of the whole stream and a calibrated error; and
// merging many independent runs, the shape of every MPI production run, must
// not shorten the bins so far that the error is under-reported.
#include "streams.hpp"

#include <alps/alea.h>
#include <alps/alea.hpp>
#include <alps/accumulators.hpp>

#include <gtest/gtest.h>
#include <gtest/gtest-spi.h>

#include <cmath>
#include <vector>

namespace {

using statistics_test::ar1_exact_error;
using statistics_test::ar1_series;

const double phi = 0.5;
const std::size_t n = std::size_t(1) << 15;

struct estimate {
    std::size_t count;
    double mean, error;
};

// The split point: half, or a place that is not a bin boundary at any level.
class SplitMerge : public ::testing::TestWithParam<std::size_t> {
  protected:
    SplitMerge() : values(ar1_series(77, phi, n)) {
        for (double value : values)
            mean += value / n;
    }
    void expect_whole(const estimate & merged, const char * library) const {
        EXPECT_EQ(merged.count, n) << library;
        EXPECT_NEAR(merged.mean, mean, 1e-13) << library;
        // One estimate from at least 128 bins scatters by 6%; the finite bin
        // length costs less than 1% at this phi.
        EXPECT_NEAR(1.0, merged.error / ar1_exact_error(phi, n), 0.3) << library;
    }
    std::vector<double> values;
    double mean = 0;
};

}

TEST_P(SplitMerge, Alea)
{
    alps::ObservableSet first, second;
    first << alps::RealObservable("x");
    second << alps::RealObservable("x");
    for (std::size_t t = 0; t < n; ++t)
        (t < GetParam() ? first : second)["x"] << values[t];
    first << second;
    const alps::RealObsevaluator merged = first["x"];
    if (GetParam() * 2 == n) {
        expect_whole({merged.count(), merged.mean(), merged.error()}, "ALEA");
    } else {
        // ALEA merges runs from their full bins: the samples of an unfinished
        // bin (16 of 10000 and 240 of 22768 here) do not enter the merged
        // mean, so it is the mean of the whole stream only to a fraction of
        // the error.
        const double exact = ar1_exact_error(phi, n);
        EXPECT_NEAR(double(merged.count()), double(n), 0.01 * n);
        EXPECT_NEAR(merged.mean(), mean, 0.1 * exact);
        EXPECT_NEAR(1.0, merged.error() / exact, 0.3);
    }
}

TEST_P(SplitMerge, Accumulators)
{
    namespace acc = alps::accumulators;
    acc::accumulator_set first, second;
    first << acc::LogBinningAccumulator<double>("x");
    second << acc::LogBinningAccumulator<double>("x");
    for (std::size_t t = 0; t < n; ++t)
        (t < GetParam() ? first : second)["x"] << values[t];
    first.merge(second);
    const acc::result_set results(first);
    expect_whole({results["x"].count(), results["x"].mean<double>(), results["x"].error<double>()}, "accumulators");
}

TEST_P(SplitMerge, AleaV2Batch)
{
    alps::alea::batch_acc<double> first(1), second(1);
    for (std::size_t t = 0; t < n; ++t)
        (t < GetParam() ? first : second) << values[t];
    first << second.finalize();
    const alps::alea::batch_result<double> merged = first.finalize();
    expect_whole({merged.count(), merged.mean()[0], merged.stderror()[0]}, "ALEA v2 batch");
}

TEST_P(SplitMerge, AleaV2Autocorr)
{
    using alps::alea::autocorr_acc;
    using alps::alea::autocorr_result;
    autocorr_acc<double> first(1), second(1), whole(1);
    for (std::size_t t = 0; t < n; ++t) {
        (t < GetParam() ? first : second) << values[t];
        whole << values[t];
    }
    const autocorr_result<double> second_result = second.finalize();
    autocorr_result<double> merged;
    const auto merge = [&] {
        first << second_result;
        merged = first.finalize();
    };
    if (GetParam() * 2 < n) {
        // Known gap, to be closed with the engine repairs: an accumulator
        // that takes a result with more levels than it has itself is left
        // in a state it cannot be finalized from.
        EXPECT_NONFATAL_FAILURE(EXPECT_NO_THROW(merge()), "throw");
        return;
    }
    ASSERT_NO_THROW(merge());
    expect_whole({merged.count(), merged.mean()[0], merged.stderror()[0]}, "ALEA v2 autocorr");

    if (GetParam() * 2 == n) {
        // Split at a bin boundary of every level below the top: the per-level
        // sums of the merged accumulators are those of the whole stream.
        const autocorr_result<double> reference = whole.finalize();
        for (std::size_t level = 0; level < 14; ++level)
            EXPECT_NEAR(reference.level(level).stderror()[0], merged.level(level).stderror()[0],
                        1e-10 * reference.level(level).stderror()[0]) << "level " << level;
    }
}

INSTANTIATE_TEST_SUITE_P(Splits, SplitMerge, ::testing::Values(n / 2, std::size_t(10000), std::size_t(22768)));

namespace {

// Eight independent runs of 2^13 samples at phi = 0.9, merged; averaged over
// replicas to make the comparison sharp. The exact error of the pooled mean is
// that of one run divided by sqrt(8).
const double slow = 0.9;
const int runs = 8;
const std::size_t run_length = std::size_t(1) << 13;
const int replicas = 40;
const double pooled_error = ar1_exact_error(slow, run_length) / std::sqrt(double(runs));
const double floor_ = 0.9;

std::uint64_t run_seed(int replica, int run) { return 104729u * replica + 31u * run + 7; }

}

TEST(ManyRuns, Alea)
{
    double error = 0;
    for (int r = 0; r < replicas; ++r) {
        alps::ObservableSet pooled;
        for (int k = 0; k < runs; ++k) {
            alps::ObservableSet run;
            run << alps::RealObservable("x");
            for (double value : ar1_series(run_seed(r, k), slow, run_length))
                run["x"] << value;
            pooled << run;
        }
        error += alps::RealObsevaluator(pooled["x"]).error() / replicas;
    }
    EXPECT_GE(error / pooled_error, floor_);
    EXPECT_LE(error / pooled_error, 1.05);
}

TEST(ManyRuns, Accumulators)
{
    namespace acc = alps::accumulators;
    double error = 0;
    for (int r = 0; r < replicas; ++r) {
        acc::accumulator_set pooled;
        pooled << acc::LogBinningAccumulator<double>("x");
        for (int k = 0; k < runs; ++k) {
            acc::accumulator_set run;
            run << acc::LogBinningAccumulator<double>("x");
            for (double value : ar1_series(run_seed(r, k), slow, run_length))
                run["x"] << value;
            pooled.merge(run);
        }
        error += acc::result_set(pooled)["x"].error<double>() / replicas;
    }
    EXPECT_GE(error / pooled_error, floor_);
    EXPECT_LE(error / pooled_error, 1.05);
}

TEST(ManyRuns, AleaV2)
{
    double autocorr_error = 0, batch_error = 0;
    for (int r = 0; r < replicas; ++r) {
        alps::alea::autocorr_acc<double> pooled_autocorr(1);
        alps::alea::batch_acc<double> pooled_batch(1);
        for (int k = 0; k < runs; ++k) {
            alps::alea::autocorr_acc<double> autocorr(1);
            alps::alea::batch_acc<double> batch(1);
            for (double value : ar1_series(run_seed(r, k), slow, run_length)) {
                autocorr << value;
                batch << value;
            }
            pooled_autocorr << autocorr.finalize();
            pooled_batch << batch.finalize();
        }
        autocorr_error += pooled_autocorr.finalize().stderror()[0] / replicas;
        batch_error += pooled_batch.finalize().stderror()[0] / replicas;
    }
    EXPECT_GE(autocorr_error / pooled_error, floor_);
    EXPECT_LE(autocorr_error / pooled_error, 1.05);
    // Known gap, to be closed with the engine repairs: merged batches keep the
    // length they had in one run, 32 samples, however many runs are merged.
    EXPECT_NONFATAL_FAILURE(EXPECT_GE(batch_error / pooled_error, floor_), "floor_");
}
