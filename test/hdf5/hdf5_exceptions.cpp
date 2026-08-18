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

#include <boost/random.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

int main() {
    std::string const filename = "test_hdf5_exceptions.h5";
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
        } catch (std::exception& ex) {
            std::string str = ex.what();
            std::size_t start = str.find_first_of("\n");
            std::cout << str.substr(0, start) << std::endl;
        }
    }
    boost::filesystem::remove(boost::filesystem::path(filename));
    return 0;
}
