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
    ASSERT_EQ(reduced.nlevel(), expected.nlevel());
    for (std::size_t level = 0; level + 1 < expected.nlevel(); ++level)
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
