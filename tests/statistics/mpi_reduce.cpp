// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Reduction over MPI ranks against a serial merge of the same slices: the
// reduction only changes how the partial sums travel, so the result on the
// root must be the serial one to rounding.
#include "streams.hpp"

#include <alps/alea.hpp>
#include <alps/alea/observables.hpp>
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

TEST(MpiReduce, AleaAutocorr)
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

TEST(MpiReduce, AleaBatch)
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
    // The i-th batches of the ranks add up, where a serial merge appends them:
    // the two estimates are equally good, not equal.
    EXPECT_NEAR(1.0, reduced.stderror()[0] / expected.stderror()[0], 0.2);
    expect_calibrated(reduced.stderror()[0], s);
    // The reduction keeps the batches of one run, at any number of ranks.
    const batch_result<double> one = fill<batch_acc<double> >(s, 0).result();
    EXPECT_EQ(reduced.store().num_batches(), one.store().num_batches());
    EXPECT_EQ(reduced.store().batch_size(), one.store().batch_size());
}

// Ranks with runs of different lengths, as after a time limit: rank r holds
// 2^(15 - 2r) samples, so the ranks differ in their number of levels and in
// their batch sizes. No sample may be lost, at any level.
TEST(MpiReduce, AleaUnequalRuns)
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
    // The longest run alone fills all 256 batches with 128 samples each. The
    // shorter runs add to its first batches; none is coarsened.
    EXPECT_EQ(batch_reduced.store().num_batches(), 256u);
    EXPECT_EQ(batch_reduced.store().batch_size(), 128u);
    EXPECT_EQ(batch_reduced.store().count()(0), 128u * s.size);
    EXPECT_EQ(batch_reduced.store().count()(255), 128u);
}

// The observable layer: a set reduced over the ranks is the set of the runs
// merged one after the other. Rank 0 has not measured "late" at all.
TEST(MpiReduce, ObservableSet)
{
    using namespace alps::alea;
    const slices s;
    const auto fill_set = [&s](int rank) {
        observable_set set;
        set << observable("x") << observable("plain").without_bins() << observable("late") << histogram("h", 0, 1);
        for (std::size_t t = 0; t < s.length; ++t) {
            const double value = s.begin(rank)[t];
            set["x"] << value;
            set["plain"] << value;
            set["h"] << (value > 0 ? 1 : 0);
            if (rank != 0)
                set["late"] << std::vector<double>({value, -value});
        }
        return set;
    };
    result_set reduced = fill_set(s.comm.rank()).results();
    reduced.reduce(mpi_reducer(s.comm, 0));
    if (s.comm.rank() != 0) {
        EXPECT_TRUE(reduced.names().empty());
        return;
    }

    observable_set serial = fill_set(0);
    for (std::size_t rank = 1; rank < s.size; ++rank)
        serial.merge(fill_set(rank).results());
    const result_set expected = serial.results();

    EXPECT_EQ(reduced["x"].count(), s.length * s.size);
    EXPECT_NEAR(reduced["x"].mean(), s.mean, 1e-13);
    EXPECT_NEAR(1.0, reduced["x"].error() / expected["x"].error(), 0.2);
    expect_calibrated(reduced["x"].error(), s);
    EXPECT_NEAR(reduced["plain"].error(), expected["plain"].error(), 1e-10 * expected["plain"].error());
    EXPECT_EQ(reduced["h"].histogram(), expected["h"].histogram());
    EXPECT_EQ(reduced["late"].count(), s.length * (s.size - 1));
    EXPECT_NEAR(reduced["late"].mean_vector()[1], expected["late"].mean_vector()[1], 1e-13);
    // The bins of the reduced set line up for a function of results.
    const observable_result square = reduced["x"] * reduced["x"], expected_square = expected["x"] * expected["x"];
    EXPECT_NEAR(1.0, square.error() / expected_square.error(), 0.2);
}

// A set's results are reduced in place unless something else holds them; then
// they are copied and stay as they were. Both ways give the same numbers.
TEST(MpiReduce, SharedResultsAreCopied)
{
    using namespace alps::alea;
    const slices s;
    observable_set set;
    set << observable("x") << observable("plain").without_bins();
    for (std::size_t t = 0; t < s.length; ++t) {
        set["x"] << s.begin(s.comm.rank())[t];
        set["plain"] << std::vector<double>({s.begin(s.comm.rank())[t], 1.0});
    }
    const result_set before = set.results();
    result_set shared = before;
    shared.reduce(mpi_reducer(s.comm, 0));
    result_set alone = set.results();
    alone.reduce(mpi_reducer(s.comm, 0));

    const result_set fresh = set.results();
    for (const char *name : {"x", "plain"}) {
        EXPECT_EQ(before[name].count(), s.length);
        EXPECT_EQ(before[name].mean_vector(), fresh[name].mean_vector());
        EXPECT_EQ(before[name].error_vector(), fresh[name].error_vector());
    }
    if (s.comm.rank() != 0)
        return;
    for (const char *name : {"x", "plain"}) {
        EXPECT_EQ(alone[name].count(), s.length * s.size);
        EXPECT_EQ(alone[name].mean_vector(), shared[name].mean_vector());
        EXPECT_EQ(alone[name].error_vector(), shared[name].error_vector());
    }
}

// The reducer sends a large observable in pieces of 4 MiB; every piece lands
// in its place.
TEST(MpiReduce, LargeObservableInPieces)
{
    using namespace alps::alea;
    alps::mpi::communicator comm;
    const std::size_t n = (std::size_t(1) << 19) + 5;   // just over one piece of doubles
    observable_set set;
    set << observable("big").without_bins();
    std::vector<double> x(n);
    for (int sample = 0; sample < 2; ++sample) {
        for (std::size_t i = 0; i < n; ++i)
            x[i] = comm.rank() + 1e-6 * i + sample;
        set["big"] << x;
    }
    result_set results = set.results();
    results.reduce(mpi_reducer(comm, 0));
    if (comm.rank() != 0)
        return;
    const std::vector<double> mean = results["big"].mean_vector();
    ASSERT_EQ(mean.size(), n);
    for (std::size_t i : {std::size_t(0), (std::size_t(1) << 19) - 1, std::size_t(1) << 19, n - 1})
        EXPECT_NEAR(mean[i], 0.5 * (comm.size() - 1) + 1e-6 * i + 0.5, 1e-12) << i;
}
