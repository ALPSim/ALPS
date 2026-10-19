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
#include <alps/parameter_input.hpp>

#include <boost/filesystem.hpp>

#include <exception>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <output.h5> < <parameter file>" << std::endl;
        return 1;
    }
    std::string const text((std::istreambuf_iterator<char>(std::cin)), std::istreambuf_iterator<char>());
    alps::Parameters parms;
    try {
        std::istringstream in(text);
        parms = alps::read_parameter_text(in);
    } catch (std::exception const & error) {
        // The text holds { } blocks when it reads as a list of runs.
        std::size_t runs = 0;
        try {
            std::istringstream in(text);
            runs = alps::read_parameter_list_text(in).size();
        } catch (std::exception const &) {}
        if (runs > 0)
            std::cerr << "p2h5 writes one parameter set, but the input has { } blocks ("
                      << runs << (runs == 1 ? " run" : " runs") << "). Use parameter2hdf5 for it." << std::endl;
        else
            std::cerr << "p2h5: " << error.what() << std::endl;
        return 1;
    }
    if (boost::filesystem::exists(boost::filesystem::path(argv[1])))
        boost::filesystem::remove(boost::filesystem::path(argv[1]));
    alps::hdf5::archive ar(argv[1], "w");
    ar["/parameters"] << parms;
}
