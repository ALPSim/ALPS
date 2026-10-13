/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Lukas Gamper <gamperl -at- gmail.com>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: abstract_task.C 3822 2010-01-30 22:02:39Z troyer $ */

#include <alps/hdf5/archive.hpp>
#include <alps/parameter.h>

#include <boost/filesystem.hpp>

#include <iostream>

int main(int argc, char **argv) {
    if (argc < 2)
        throw std::invalid_argument("no name passed");
    alps::Parameters parms;
    std::cin >> parms;
    if (boost::filesystem::exists(boost::filesystem::path(argv[1])))
        boost::filesystem::remove(boost::filesystem::path(argv[1]));
    alps::hdf5::archive ar(argv[1], "w");
    ar["/parameters"] << parms;
}
