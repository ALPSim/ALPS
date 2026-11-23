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
#include <thread>
#include <vector>
#include <iostream>
#include <algorithm>

TEST(Hdf5, Multiarchive) {
    alps::testing::TemporaryDirectory temporary;
    std::string const filename = (temporary.path() / "test_hdf5_multiarchive.h5").string();
    if (std::filesystem::exists(std::filesystem::path(filename)))
        std::filesystem::remove(std::filesystem::path(filename));
    {
        using namespace alps;
        alps::hdf5::archive oar(filename, "a");
    }
    {
        using namespace alps;
        alps::hdf5::archive iar(filename);
        alps::hdf5::archive oar(filename, "a");
        oar << make_pvp("/data", 42);
    }
    {
        using namespace alps;
        alps::hdf5::archive iar(filename, "r");
        int test;
        iar >> make_pvp("/data", test);
        EXPECT_EQ(test, 42);
        {
            alps::hdf5::archive iar2(filename, "r");
            int test2;
            iar2 >> make_pvp("/data", test2);
            EXPECT_EQ(test2, 42);
            iar >> make_pvp("/data", test);
            EXPECT_EQ(test, 42);
        }
        iar >> make_pvp("/data", test);
        EXPECT_EQ(test, 42);
        {
            alps::hdf5::archive iar3(filename, "r");
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
        alps::hdf5::archive iar4(filename, "r");
        int test4;
        iar4 >> make_pvp("/data", test4);
        EXPECT_EQ(test4, 42);
    }
    std::filesystem::remove(std::filesystem::path(filename));

}

// One thread reads a file while another opens it for writing. Opening it for
// writing reopens the shared file; with a thread-safe HDF5 that used to abort
// with "Not all resources closed".
TEST(Hdf5, ReaderAndWriterInTwoThreads) {
    alps::testing::TemporaryDirectory temporary;
    std::string const filename = (temporary.path() / "mixed.h5").string();
    { alps::hdf5::archive ar(filename, "w"); ar["/x"] << 1; }
    std::vector<int> read(500, 0);
    std::thread reader([&] {
        for (int& x : read) { alps::hdf5::archive ar(filename, "r"); ar["/x"] >> x; }
    });
    std::thread writer([&] {
        for (int i = 0; i < 500; ++i) { alps::hdf5::archive ar(filename, "a"); ar["/y" + std::to_string(i)] << i; }
    });
    reader.join();
    writer.join();
    EXPECT_EQ(std::count(read.begin(), read.end(), 1), 500);
    alps::hdf5::archive ar(filename, "r");
    int last = -1;
    ar["/y499"] >> last;
    EXPECT_EQ(last, 499);
}
