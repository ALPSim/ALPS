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
#include <iostream>
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/multi_array.hpp>

#include <boost/multi_array.hpp>
#include <boost/filesystem.hpp>

#include <vector>

using namespace std;
using boost::multi_array;

TEST(Hdf5, MultiArray)
{
    alps::testing::TemporaryDirectory temporary;
    if (boost::filesystem::exists(boost::filesystem::path((temporary.path() / "test_hdf5_multi_array.h5").string())))
        boost::filesystem::remove(boost::filesystem::path((temporary.path() / "test_hdf5_multi_array.h5").string()));

    multi_array<double,2> a( boost::extents[3][3] );
    multi_array<double,2> b( boost::extents[4][4] );

    std::fill(a.data(), a.data() + a.num_elements(), 1.25);
    std::fill(b.data(), b.data() + b.num_elements(), -2.5);
    // Write
    {
        alps::hdf5::archive ar((temporary.path() / "test_hdf5_multi_array.h5").string(),"a");
        vector< multi_array<double,2> > v(2,a);
        ar << alps::make_pvp("uniform",v);
        v.push_back(b);
        ar << alps::make_pvp("nonuniform",v);
    }

    // Read
    {
        alps::hdf5::archive ar((temporary.path() / "test_hdf5_multi_array.h5").string(),"r");
        vector< multi_array<double,2> > w;
        ar >> alps::make_pvp("nonuniform",w);
        ASSERT_EQ(w.size(), 3u);
        EXPECT_EQ(w[0], a);
        EXPECT_EQ(w[1], a);
        EXPECT_EQ(w[2], b);
        ar >> alps::make_pvp("uniform",w);
        ASSERT_EQ(w.size(), 2u);
        EXPECT_EQ(w[0], a);
        EXPECT_EQ(w[1], a);
    }

    if (boost::filesystem::exists(boost::filesystem::path((temporary.path() / "test_hdf5_multi_array.h5").string())))
        boost::filesystem::remove(boost::filesystem::path((temporary.path() / "test_hdf5_multi_array.h5").string()));

}
