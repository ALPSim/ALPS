/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "sum_simulation.hpp"
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <alps/check_schedule.hpp>
#include <new>

TEST(MonteCarloSchedule, ChecksImmediatelyAfterConstructionInReusedStorage) {
    // A stale NaN interval prevented MPI runs from ever checking completion.
    alignas(alps::check_schedule) unsigned char storage[sizeof(alps::check_schedule)];
    volatile unsigned char* bytes = storage;
    for (std::size_t i = 0; i < sizeof(storage); ++i) bytes[i] = 0xff;
    auto* schedule = new (storage) alps::check_schedule;
    EXPECT_DOUBLE_EQ(schedule->check_interval(), 0.);
    EXPECT_TRUE(schedule->pending());
    schedule->update(0.25);
    EXPECT_DOUBLE_EQ(schedule->check_interval(), 60.);
    EXPECT_FALSE(schedule->pending());
    schedule->~check_schedule();
}

TEST(MonteCarloRunner, CompletesAndPersistsScalarAndVectorResults) {
    alps::params parameters;
    parameters["COUNT"] = 256;
    parameters["SEED"] = 42;
    sum_simulation simulation(parameters);
    EXPECT_DOUBLE_EQ(simulation.fraction_completed(), 0.);
    ASSERT_TRUE(simulation.run([] { return false; }));
    EXPECT_EQ(simulation.count, 256);
    EXPECT_DOUBLE_EQ(simulation.fraction_completed(), 1.);
    const auto results = alps::collect_results(simulation);
    const double expected = simulation.sum / simulation.count;
    EXPECT_EQ(results["SValue"].count(), 256u);
    // Only floating-point accumulation order differs from the explicit sum.
    EXPECT_NEAR(results["SValue"].mean<double>(), expected, 1e-12);
    EXPECT_TRUE(std::isfinite(results["SValue"].error<double>()));
    const auto vector = results["VValue"].mean<std::vector<double>>();
    ASSERT_EQ(vector.size(), 3u);
    for (double value : vector) EXPECT_NEAR(value, expected, 1e-12);
    EXPECT_GE(expected, std::exp(-1.));
    EXPECT_LE(expected, 1.);

    alps::testing::TemporaryDirectory directory;
    const auto filename = (directory.path() / "results.h5").string();
    alps::save_results(results, parameters, filename, "/simulation/results");
    alps::hdf5::archive archive(filename);
    alps::mcresults restored;
    archive["/simulation/results"] >> restored;
    EXPECT_EQ(restored["SValue"].count(), results["SValue"].count());
    EXPECT_DOUBLE_EQ(restored["SValue"].mean<double>(), results["SValue"].mean<double>());
    EXPECT_EQ(restored["VValue"].mean<std::vector<double>>(), vector);
}

TEST(MonteCarloRunner, StopCallbackPreventsFurtherMeasurements) {
    alps::params parameters;
    parameters["COUNT"] = 256;
    parameters["SEED"] = 42;
    sum_simulation simulation(parameters);
    EXPECT_FALSE(simulation.run([&] { return simulation.count == 32; }));
    EXPECT_EQ(simulation.count, 32);
    EXPECT_DOUBLE_EQ(simulation.fraction_completed(), 0.125);
    EXPECT_EQ(alps::collect_results(simulation)["SValue"].count(), 32u);
}
