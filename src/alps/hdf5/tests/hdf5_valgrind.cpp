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
#include <sstream>
#include <vector>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>

TEST(Hdf5, Valgrind) {
    alps::testing::TemporaryDirectory temporary;

    for (int i=0; i<100; ++i) {
        std::vector<double> vec(10,2.);
        alps::hdf5::archive ar((temporary.path() / "test_hdf5_valgrind.h5").string(), "w");
        std::ostringstream ss;
        ss << "/vec" << i;
        SCOPED_TRACE(i);
        ar << alps::make_pvp(ss.str(), vec);
        std::vector<double> restored;
        ar[ss.str()] >> restored;
        EXPECT_EQ(restored, vec);
    }

}
