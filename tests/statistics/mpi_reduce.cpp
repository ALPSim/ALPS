// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Reduction over MPI ranks against a serial merge of the same slices: the
// reduction only changes how the partial sums travel, so the result on the
// root must be the serial one to rounding.
#include "streams.hpp"

#include <alps/alea.hpp>
#include <alps/accumulators.hpp>
#include <alps/utilities/mpi.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>
#include <vector>

namespace {

const double phi = 0.5;
const std::size_t total = 3 * (std::size_t(1) << 14);

// Rank r of size ranks measures the r-th contiguous slice of one stream.
struct slices {
    slices() : values(statistics_test::ar1_series(4711, phi, total)), size(comm.size()), length(total / size) {
        for (std::size_t t = 0; t < length * size; ++t)
            mean += values[t] / (length * size);
    }
    const double * begin(int rank) const { return values.data() + rank * length; }
    alps::mpi::communicator comm;
    std::vector<double> values;
    std::size_t size, length;
    double mean = 0;
};

template <typename Acc> Acc fill(const slices & s, int rank) {
    Acc acc(1);
    for (std::size_t t = 0; t < s.length; ++t)
        acc << s.begin(rank)[t];
    return acc;
}

void expect_calibrated(double error, const slices & s) {
    // One estimate from at least 128 bins, as in merge.cpp.
    EXPECT_NEAR(1.0, error / statistics_test::ar1_exact_error(phi, s.length * s.size), 0.3);
}

}

TEST(MpiReduce, AleaV2Autocorr)
{
    using namespace alps::alea;
    const slices s;
    autocorr_result<double> reduced = fill<autocorr_acc<double> >(s, s.comm.rank()).result();
    reduced.reduce(mpi_reducer(s.comm, 0));
    if (s.comm.rank() != 0)
        return;

    autocorr_acc<double> serial = fill<autocorr_acc<double> >(s, 0);
    for (std::size_t rank = 1; rank < s.size; ++rank)
        serial << fill<autocorr_acc<double> >(s, rank).result();
    const autocorr_result<double> expected = serial.finalize();

    EXPECT_EQ(reduced.count(), s.length * s.size);
    EXPECT_NEAR(reduced.mean()[0], s.mean, 1e-13);
    // Merging one after the other adds levels as the count grows; the
    // reduction keeps the levels of the longest run.
    ASSERT_LE(reduced.nlevel(), expected.nlevel());
    for (std::size_t level = 0; level + 1 < reduced.nlevel(); ++level)
        EXPECT_NEAR(expected.level(level).stderror()[0], reduced.level(level).stderror()[0],
                    1e-10 * expected.level(level).stderror()[0]) << "level " << level;
    expect_calibrated(reduced.stderror()[0], s);
}

TEST(MpiReduce, AleaV2Batch)
{
    using namespace alps::alea;
    const slices s;
    batch_result<double> reduced = fill<batch_acc<double> >(s, s.comm.rank()).result();
    reduced.reduce(mpi_reducer(s.comm, 0));
    if (s.comm.rank() != 0)
        return;

    batch_acc<double> serial = fill<batch_acc<double> >(s, 0);
    for (std::size_t rank = 1; rank < s.size; ++rank)
        serial << fill<batch_acc<double> >(s, rank).result();
    const batch_result<double> expected = serial.finalize();

    EXPECT_EQ(reduced.count(), s.length * s.size);
    EXPECT_NEAR(reduced.mean()[0], s.mean, 1e-13);
    EXPECT_NEAR(expected.stderror()[0], reduced.stderror()[0], 1e-10 * expected.stderror()[0]);
    expect_calibrated(reduced.stderror()[0], s);
}

// Ranks with runs of different lengths, as after a time limit: rank r holds
// 2^(15 - 2r) samples, so the ranks differ in their number of levels and in
// their batch sizes. No sample may be lost, at any level.
TEST(MpiReduce, AleaV2UnequalRuns)
{
    using namespace alps::alea;
    const slices s;
    const auto length = [](std::size_t rank) { return std::size_t(1) << (15 - 2 * rank); };
    std::size_t total_length = 0;
    double sum = 0;
    for (std::size_t rank = 0; rank < s.size; ++rank) {
        total_length += length(rank);
        sum = std::accumulate(s.begin(rank), s.begin(rank) + length(rank), sum);
    }
    autocorr_acc<double> autocorr(1);
    batch_acc<double> batch(1);
    for (std::size_t t = 0; t < length(s.comm.rank()); ++t) {
        autocorr << s.begin(s.comm.rank())[t];
        batch << s.begin(s.comm.rank())[t];
    }
    autocorr_result<double> autocorr_reduced = autocorr.finalize();
    batch_result<double> batch_reduced = batch.finalize();
    autocorr_reduced.reduce(mpi_reducer(s.comm, 0));
    batch_reduced.reduce(mpi_reducer(s.comm, 0));
    if (s.comm.rank() != 0)
        return;

    EXPECT_NEAR(autocorr_reduced.mean()[0], sum / total_length, 1e-13);
    for (std::size_t level = 0; level < autocorr_reduced.nlevel(); ++level)
        EXPECT_EQ(autocorr_reduced.level(level).count(), total_length) << "level " << level;
    EXPECT_EQ(batch_reduced.count(), total_length);
    EXPECT_NEAR(batch_reduced.mean()[0], sum / total_length, 1e-13);
    // The longest run alone fills all 256 batches with 128 samples each. With
    // the others appended, neighbours are added once.
    EXPECT_EQ(batch_reduced.store().batch_size(), 256u);
    EXPECT_EQ(batch_reduced.store().count()(0), 256u);
}

TEST(MpiReduce, Accumulators)
{
    namespace acc = alps::accumulators;
    const slices s;
    const auto fill_set = [&s](int rank) {
        acc::accumulator_set set;
        set << acc::LogBinningAccumulator<double>("x");
        for (std::size_t t = 0; t < s.length; ++t)
            set["x"] << s.begin(rank)[t];
        return set;
    };
    acc::accumulator_set reduced = fill_set(s.comm.rank());
    reduced["x"].collective_merge(s.comm, 0);
    if (s.comm.rank() != 0)
        return;

    acc::accumulator_set serial = fill_set(0);
    for (std::size_t rank = 1; rank < s.size; ++rank)
        serial.merge(fill_set(rank));
    const acc::result_set expected(serial), results(reduced);

    EXPECT_EQ(results["x"].count(), s.length * s.size);
    EXPECT_NEAR(results["x"].mean<double>(), s.mean, 1e-13);
    EXPECT_NEAR(expected["x"].error<double>(), results["x"].error<double>(), 1e-10 * expected["x"].error<double>());
    expect_calibrated(results["x"].error<double>(), s);
}
