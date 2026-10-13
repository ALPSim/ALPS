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

#include <algorithm>
#include <vector>
#include <complex>

#include <alps/hdf5.hpp>
#include <alps/hdf5/vector.hpp>

TEST(Hdf5, Pair) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive ar((temporary.path() / "creal.h5").string(), "a");
    {
        std::vector<double> a(1e6);
        ar << alps::make_pvp("a",
            std::make_pair(
                static_cast< double const *>(&a.front())
                , std::vector<std::size_t>(1,a.size())
            )
        );
    }
    {
        std::vector<std::complex<double> > a(1e6);
        ar << alps::make_pvp("a",
            std::make_pair(
                static_cast<std::complex<double> const *>(&a.front())
                , std::vector<std::size_t>(1,a.size())
            )
        );
    }
    std::vector<std::complex<double>> restored;
    ar["a"] >> restored;
    ASSERT_EQ(restored.size(), 1000000u);
    EXPECT_TRUE(std::all_of(restored.begin(), restored.end(), [](std::complex<double> value) {
        return value == std::complex<double>();
    }));
}
