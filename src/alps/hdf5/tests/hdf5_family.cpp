/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include <alps/hdf5/archive.hpp>

#include <filesystem>
#include <boost/random.hpp>

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

TEST(Hdf5, Family) {
    alps::testing::TemporaryDirectory temporary;
    std::string const filename = (temporary.path() / "test%05d.h5").string();
    if (std::filesystem::exists(std::filesystem::path(filename)))
        std::filesystem::remove(std::filesystem::path(filename));
    {
        using namespace alps;
        alps::hdf5::archive oar(filename, "al");
    }
    {
        using namespace alps;
        alps::hdf5::archive oar(filename, "al");
        oar << make_pvp("/data", 42);
    }
    {
        using namespace alps;
        alps::hdf5::archive iar(filename, "l");
        int test;
        iar >> make_pvp("/data", test);
        EXPECT_EQ(test, 42);
        {
            alps::hdf5::archive iar2(filename, "l");
            int test2;
            iar2 >> make_pvp("/data", test2);
            EXPECT_EQ(test2, 42);
            iar >> make_pvp("/data", test);
            EXPECT_EQ(test, 42);
        }
        iar >> make_pvp("/data", test);
        EXPECT_EQ(test, 42);
        {
            alps::hdf5::archive iar3(filename, "l");
            int test3;
            iar >> make_pvp("/data", test);
            EXPECT_EQ(test, 42);
            iar3 >> make_pvp("/data", test3);
            EXPECT_EQ(test3, 42);
        }
        iar >> make_pvp("/data", test);
        EXPECT_EQ(test, 42);
    }
    {
        using namespace alps;
        alps::hdf5::archive iar4(filename, "l");
        int test4;
        iar4 >> make_pvp("/data", test4);
        EXPECT_EQ(test4, 42);
    }
    std::filesystem::remove(std::filesystem::path(filename));

}
