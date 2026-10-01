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

#include <iostream>
#include <sstream>
#include <vector>
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>

int main() {

    for (int i=0; i<100; ++i) {
        std::vector<double> vec(10,2.);
        alps::hdf5::archive ar("test_hdf5_valgrind.h5", "w");
        std::ostringstream ss;
        ss << "/vec" << i;
        ar << alps::make_pvp(ss.str(), vec);
    }
    return 0;
}
