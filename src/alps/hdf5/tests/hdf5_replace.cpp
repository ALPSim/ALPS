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

#include <vector>

using namespace std;

TEST(Hdf5, Replace) {
    alps::testing::TemporaryDirectory temporary;

   vector<double> vec(100, 10.);

   {
       alps::hdf5::archive ar((temporary.path() / "res.h5").string(), "a");
       ar << alps::make_pvp("/vec2", vec);
   }
   {
       alps::hdf5::archive ar((temporary.path() / "res.h5").string(), "w");
       ar << alps::make_pvp("/vec", vec);
   }
   {
       vector<double> tmp;
       alps::hdf5::archive ar((temporary.path() / "res.h5").string());
       ar >> alps::make_pvp("/vec2", tmp);
       EXPECT_EQ(tmp, vec);
       ar >> alps::make_pvp("/vec", tmp);
       EXPECT_EQ(tmp, vec);
   }


}
