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

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <vector>

// Explicit opt-in: the largest allocation is 2 GiB and the family uses about
// 4 GiB of disk. Ordinary extensive type tests do not enable this stress test.
TEST(Hdf5LargeArchive, FamilyRetainsLargeDatasets) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive archive((temporary.path() / "large%d.h5").string(), "al");
    for (std::size_t size = 1; size < (1ULL << 29); size <<= 1) {
        SCOPED_TRACE(size);
        const auto path = "/" + std::to_string(size);
        const std::vector<double> values(size, 10.);
        archive[path] << values;
        EXPECT_EQ(archive.extent(path), std::vector<std::size_t>({size}));
        double first = 0., last = 0.;
        archive.read(path, &first, std::vector<std::size_t>{1}, std::vector<std::size_t>{0});
        archive.read(path, &last, std::vector<std::size_t>{1}, std::vector<std::size_t>{size - 1});
        EXPECT_EQ(first, 10.);
        EXPECT_EQ(last, 10.);
    }
}
