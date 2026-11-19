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


#include <alps/alea/observables.hpp>
#include <alps/alea/observables_io.hpp>
#include <alps/hdf5.hpp>
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

TEST(Hdf5ObservableSet, ArchiveRetainsCountAndMean) {
    alps::testing::TemporaryDirectory temporary;
    const auto filename = (temporary.path() / "observables.h5").string();
    alps::alea::observable_set measurements;
    measurements << alps::alea::observable("E");
    measurements["E"] << 1.;
    {
        alps::hdf5::archive archive(filename, "w");
        alps::alea::save_results(archive, "/simulation/results", measurements.results());
    }
    alps::hdf5::archive archive(filename, "r");
    unsigned long long count = 0;
    double mean = 0.;
    archive["/simulation/results/E/count"] >> count;
    archive["/simulation/results/E/mean/value"] >> mean;
    EXPECT_EQ(count, 1u);
    EXPECT_EQ(mean, 1.);
    const alps::alea::result_set restored = alps::alea::load_results(archive, "/simulation/results");
    EXPECT_EQ(restored["E"].count(), 1u);
    EXPECT_EQ(restored["E"].mean(), 1.);
}
