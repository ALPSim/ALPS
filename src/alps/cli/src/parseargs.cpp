/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/parseargs.hpp>
#include <alps/cli/options.hpp>

#include <sstream>
#include <iostream>

namespace alps {

    parseargs::parseargs(int argc, char *argv[]) {
        cli::options options("Options");
        options.flag("continue,c", "load simulation from checkpoint")
            .value("timelimit,T", timelimit, 0, "time limit for the simulation")
            .value("Tmin,i", tmin, 1, "minimum time to check if simulation has finished")
            .value("Tmax,a", tmax, 600, "maximum time to check if simulation has finished")
            .value("inputfile", input_file, "input file in hdf5 or xml format")
            .value("outputfile", output_file, "", "output file in hdf5 format")
            .positional("inputfile").positional("outputfile");

        try {
            options.parse(argc, argv);

            resume = options.count("continue");
            if (output_file.empty())
                output_file = input_file.substr(0, input_file.find_last_of('.')) +  ".out.h5";
        } catch (...) {
    		std::stringstream ss;
            ss << "usage: [-T timelimit] [-i tmin] [-a tmax] [-c] inputfile [outputfile]" << std::endl
               << options << std::endl;
            std::cerr << ss.str();
            std::abort();
        }
    }

}
