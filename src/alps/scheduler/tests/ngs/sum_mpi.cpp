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
#include <alps/mcmpiadapter.hpp>
#include <alps/testing/temporary_directory.hpp>
#include <gtest/gtest.h>
#include <boost/mpi/collectives.hpp>
#include <functional>

// Make completion checks deterministic instead of depending on wall time.
struct every_step {
    bool pending() const { return true; }
    void update(double) {}
};

TEST(ParallelMonteCarloRunner, CollectsAllRanksAndPersistsResults) {
    boost::mpi::communicator world;
    ASSERT_GE(world.size(), 2);
    alps::params parameters;
    if (world.rank() == 0) {
        parameters["COUNT"] = 128 * world.size();
        parameters["SEED"] = 42;
    }
    alps::broadcast(world, parameters);
    alps::mcmpiadapter<sum_simulation, every_step> simulation(parameters, world);
    EXPECT_TRUE(simulation.run([] { return false; }));
    EXPECT_EQ(simulation.count, 128);
    EXPECT_DOUBLE_EQ(simulation.fraction_completed(), 1.);
    const double sum = boost::mpi::all_reduce(world, simulation.sum, std::plus<double>());
    const auto results = alps::collect_results(simulation);
    // Every rank must finish the collectives before rank-specific assertions.
    if (world.rank() == 0) {
        const double expected = sum / (128 * world.size());
        EXPECT_EQ(results["SValue"].count(), 128u * world.size());
        EXPECT_NEAR(results["SValue"].mean<double>(), expected, 1e-12);
        EXPECT_TRUE(std::isfinite(results["SValue"].error<double>()));
        const auto vector = results["VValue"].mean<std::vector<double>>();
        ASSERT_EQ(vector.size(), 3u);
        for (double value : vector) EXPECT_NEAR(value, expected, 1e-12);
        EXPECT_NEAR((results["SValue"] + 1.).mean<double>(), expected + 1., 1e-12);
        EXPECT_NEAR((results["SValue"] + results["SValue"]).mean<double>(), 2. * expected, 1e-12);
        alps::testing::TemporaryDirectory directory;
        const auto filename = (directory.path() / "results.h5").string();
        alps::save_results(results, parameters, filename, "/simulation/results");
        alps::hdf5::archive archive(filename);
        alps::mcresults restored;
        archive["/simulation/results"] >> restored;
        EXPECT_DOUBLE_EQ(restored["SValue"].mean<double>(), results["SValue"].mean<double>());
        EXPECT_EQ(restored["VValue"].mean<std::vector<double>>(), vector);
    }
}
