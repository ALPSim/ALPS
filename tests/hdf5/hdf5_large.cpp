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

#include <alps/hdf5/archive.hpp>

#include <boost/filesystem.hpp>

#include <vector>

using namespace std;
using namespace alps;

int main () {

    for (std::size_t i = 0; i < 100; ++i)
        if (boost::filesystem::exists(boost::filesystem::path("large" + alps::cast<std::string>(i) + ".h5")))
            boost::filesystem::remove(boost::filesystem::path("large" + alps::cast<std::string>(i) + ".h5"));

    hdf5::archive ar("large%d.h5", "al");
    for (unsigned long long s = 1; s < (1ULL << 29); s <<= 1) {
        std::cout << s << std::endl;
        vector<double> vec(s, 10.);
        ar << make_pvp("/" + cast<std::string>(s), vec);
    }

    for (std::size_t i = 0; i < 100; ++i)
        if (boost::filesystem::exists(boost::filesystem::path("large" + alps::cast<std::string>(i) + ".h5")))
            boost::filesystem::remove(boost::filesystem::path("large" + alps::cast<std::string>(i) + ".h5"));
    return 0;
}
