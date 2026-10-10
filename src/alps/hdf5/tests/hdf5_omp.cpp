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

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include <iostream>
#include <alps/hdf5.hpp>

#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

TEST(Hdf5OpenMP, ConcurrentArchivesRetainEachWorkersValue) {
    alps::testing::TemporaryDirectory temporary;
    std::vector<std::string> failures(10);
#pragma omp parallel for
    for (int i = 0; i < 10; ++i) {
        try {
            const auto filename = (temporary.path() / ("omp." + std::to_string(i) + ".h5")).string();
            { alps::hdf5::archive archive(filename, "w"); archive["/value"] << i; }
            alps::hdf5::archive archive(filename, "r");
            int restored = -1;
            archive["/value"] >> restored;
            if (restored != i) failures[i] = "Restored value differs from worker index";
        } catch (const std::exception& error) {
            failures[i] = error.what();
        }
    }
    for (int i = 0; i < 10; ++i) EXPECT_TRUE(failures[i].empty()) << "Worker " << i << ": " << failures[i];
}
