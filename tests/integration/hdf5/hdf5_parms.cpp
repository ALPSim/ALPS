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
#include <boost/filesystem.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/parameter.h>

using namespace std;

TEST(Hdf5LegacyParameters, OverwriteChangesValuesAndTypes)
{
    alps::testing::TemporaryDirectory temporary;
    if (boost::filesystem::exists((temporary.path() / "parms.h5").string()) && boost::filesystem::is_regular_file((temporary.path() / "parms.h5").string()))
        boost::filesystem::remove((temporary.path() / "parms.h5").string());

    alps::Parameters p, p2;
    p["a"] = 10;
    p["b"] = "test";
    p["c"] = 10.;
    p["d"] = 5.;

    {
        alps::hdf5::archive ar((temporary.path() / "parms.h5").string(), "a");
        ar << alps::make_pvp("/parameters", p);
    }
    {
        alps::hdf5::archive ar((temporary.path() / "parms.h5").string(), "r");
        alps::Parameters pin;
        ar >> alps::make_pvp("/parameters", pin);
        EXPECT_EQ(std::string(pin["a"]), "10");
        EXPECT_EQ(std::string(pin["b"]), "test");
        EXPECT_EQ(std::string(pin["c"]), "10");
        EXPECT_EQ(std::string(pin["d"]), "5");
    }

    // "a" is modified from int to double
    // "c" is modified from double to double (but with decimals)
    p2["a"] = 10.5;
    p2["c"] = 5.2;
    {
        alps::hdf5::archive ar((temporary.path() / "parms.h5").string(), "a");
        ar << alps::make_pvp("/parameters", p2);
    }
    {
        alps::hdf5::archive ar((temporary.path() / "parms.h5").string(), "r");
        alps::Parameters pin;
        ar >> alps::make_pvp("/parameters", pin);
        EXPECT_EQ(std::stod(std::string(pin["a"])), 10.5);
        EXPECT_EQ(std::stod(std::string(pin["c"])), 5.2);
        EXPECT_EQ(std::string(pin["b"]), "test");
        EXPECT_EQ(std::string(pin["d"]), "5");
    }

    // "d" is modified from double to string
    p2["d"] = "newtype";
    {
        alps::hdf5::archive ar((temporary.path() / "parms.h5").string(), "a");
        ar << alps::make_pvp("/parameters", p2);
    }
    {
        alps::hdf5::archive ar((temporary.path() / "parms.h5").string(), "r");
        alps::Parameters pin;
        ar >> alps::make_pvp("/parameters", pin);
        EXPECT_EQ(std::stod(std::string(pin["a"])), 10.5);
        EXPECT_EQ(std::stod(std::string(pin["c"])), 5.2);
        EXPECT_EQ(std::string(pin["d"]), "newtype");
    }


}
