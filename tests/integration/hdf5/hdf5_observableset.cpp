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


#include <alps/alea.h>
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

TEST(Hdf5ObservableSet, ArchiveRetainsCountAndMean) {
    alps::testing::TemporaryDirectory temporary;
    const auto filename = (temporary.path() / "observables.h5").string();
    alps::ObservableSet measurements;
    measurements << alps::RealObservable("E");
    measurements.get<alps::RealObservable>("E") << 1;
    {
        alps::hdf5::archive archive(filename, "w");
        archive["/simulation/results"] << measurements;
    }
    alps::hdf5::archive archive(filename, "r");
    unsigned long long count = 0;
    double mean = 0.;
    archive["/simulation/results/E/count"] >> count;
    archive["/simulation/results/E/mean/value"] >> mean;
    EXPECT_EQ(count, 1u);
    EXPECT_EQ(mean, 1.);
    alps::ObservableSet restored;
    archive["/simulation/results"] >> restored;
    const auto& energy = restored.get<alps::RealObservable>("E");
    EXPECT_EQ(energy.count(), 1u);
    EXPECT_EQ(energy.mean(), 1.);
}
