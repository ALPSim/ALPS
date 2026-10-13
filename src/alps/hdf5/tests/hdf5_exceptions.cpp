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

#include <alps/hdf5/archive.hpp>

#include <boost/random.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <algorithm>

TEST(Hdf5, Exceptions) {
    alps::testing::TemporaryDirectory temporary;
    std::string const filename = (temporary.path() / "test_hdf5_exceptions.h5").string();
    if (boost::filesystem::exists(boost::filesystem::path(filename)))
        boost::filesystem::remove(boost::filesystem::path(filename));
    {
        alps::hdf5::archive oar(filename, "a");
    }
    {
        using namespace alps;
        alps::hdf5::archive iar(filename, "r");
        double test;
        try {
            iar >> make_pvp("/not/existing/path", test);
            FAIL() << "Reading a missing dataset did not throw";
        } catch (std::exception& ex) {
            std::string str = ex.what();
            std::size_t start = str.find_first_of("\n");
            EXPECT_EQ(str.substr(0, start), "the path does not exist: /not/existing/path");
        }
    }
    boost::filesystem::remove(boost::filesystem::path(filename));

}

TEST(Hdf5, StringConversionFailureKeepsArchiveUsable) {
    alps::testing::TemporaryDirectory temporary;
    alps::hdf5::archive archive((temporary.path() / "conversion.h5").string(), "w");
    archive.create_group("/scalar");
    std::string const overflow = "999999999999999999999999999999999999999999999999";

    for (std::string const path : {"/scalar/value", "/scalar/@value"}) {
        SCOPED_TRACE(path);
        archive[path] << overflow;

        int number = 17;
        EXPECT_THROW(archive[path] >> number, std::out_of_range);
        EXPECT_EQ(number, 17);

        std::string restored;
        archive[path] >> restored;
        EXPECT_EQ(restored, overflow);
    }
}
