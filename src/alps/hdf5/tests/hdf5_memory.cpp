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
#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>

#include <boost/filesystem.hpp>

#include <iostream>
using namespace std;

TEST(Hdf5, Memory)
{
    alps::testing::TemporaryDirectory temporary;
    if (boost::filesystem::exists(boost::filesystem::path((temporary.path() / "test_hdf5_memory.h5").string())))
        boost::filesystem::remove(boost::filesystem::path((temporary.path() / "test_hdf5_memory.h5").string()));
    {
        alps::hdf5::archive oa((temporary.path() / "test_hdf5_memory.h5").string(), "w");
        std::vector<std::complex<double> > foo(3);
        std::vector<double> foo2(3);
        oa << alps::make_pvp("/foo", foo);
        oa << alps::make_pvp("/foo2", foo2);
    }

    {
        alps::hdf5::archive archive((temporary.path() / "test_hdf5_memory.h5").string());
        std::vector<double> real;
        EXPECT_THROW(archive["/foo"] >> real, alps::hdf5::archive_error);
        archive["/foo2"] >> real;
        EXPECT_EQ(real, std::vector<double>(3, 0.));
    }
}
