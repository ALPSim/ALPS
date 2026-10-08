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
#include <alps/hdf5/vector.hpp>
#include <alps/hdf5/complex.hpp>

#include <boost/filesystem.hpp>

#include <vector>
#include <complex>

using namespace std;

TEST(Hdf5, Vecveccplx)
{
    alps::testing::TemporaryDirectory temporary;
    if (boost::filesystem::exists(boost::filesystem::path((temporary.path() / "vvcplx.h5").string())))
        boost::filesystem::remove(boost::filesystem::path((temporary.path() / "vvcplx.h5").string()));
    {
      vector< vector< complex<double> > > v;
      for( int i = 0; i < 3; ++i )
        v.push_back(vector< complex<double> >(i+1, complex<double>(i,2*i)));
      alps::hdf5::archive ar((temporary.path() / "vvcplx.h5").string(),alps::hdf5::archive::WRITE);
      ar << alps::make_pvp("v",v);
      vector<vector<complex<double>>> restored;
      ar["v"] >> restored;
      EXPECT_EQ(restored, v);
    }
    boost::filesystem::remove(boost::filesystem::path((temporary.path() / "vvcplx.h5").string()));
}
