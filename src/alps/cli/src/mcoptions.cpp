/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/mcoptions.hpp>
#include <alps/cli.hpp>
#include <alps/cli/options.hpp>
#include <alps/utilities/stacktrace.hpp>

#include <iostream>
#include <stdexcept>

namespace alps {

    mcoptions::mcoptions(int argc, char* argv[]) : mcoptions(argc, argv, "framework") {}

    mcoptions::mcoptions(int argc, char* argv[], const std::string& citation_component) : valid(false), resume(false), type(SINGLE) {
        cli::options desc("Allowed options");
        desc.flag("help", "produce help message")
            .flag("single", "run single process")
            #ifndef ALPS_NGS_SINGLE_THREAD
                .flag("threaded", "run in multithread environment")
            #endif
            #ifdef ALPS_HAVE_MPI
                .flag("mpi", "run in parallel using MPI")
            #endif
            .flag("continue,c", "load simulation from checkpoint")
            .value("time-limit,T", time_limit, 0, "time limit for the simulation")
            .value("input-file", input_file, "input file in hdf5 format")
            .value("output-file", output_file, "<unspecified>", "output file in hdf5 format")
            .value("checkpoint-file", checkpoint_file, "", "checkpoint file in hdf5 format")
            .positional("input-file").positional("output-file").positional("checkpoint-file");
        if (handle_cli_information(argc, argv, citation_component, [&] { std::cout << desc << "\n"; })) return;
        desc.parse(argc, argv);
        valid = true;
        if (input_file.empty())
            throw std::invalid_argument("No job file specified" + ALPS_STACKTRACE);
        if (desc.count("threaded") && desc.count("mpi"))
            type = HYBRID;
        else if (desc.count("threaded"))
            #ifdef ALPS_NGS_SINGLE_THREAD
                throw std::logic_error("Not build with multithread support" + ALPS_STACKTRACE);
            #else
                type = THREADED;
            #endif
        else if (desc.count("mpi")) {
            type = MPI;
            #ifndef ALPS_HAVE_MPI
                throw std::logic_error("Not build with MPI" + ALPS_STACKTRACE);
            #endif
        }
        if (desc.count("continue"))
            resume = true;
        if (output_file == "<unspecified>") {
            if (input_file.find(".in.h5") != std::string::npos)
              output_file = input_file.substr(0,input_file.find_last_of(".in.h5")-5)+ ".out.h5";
            else if (input_file.find(".out.h5") != std::string::npos)
              output_file = input_file;
            else
              output_file = input_file.substr(0,input_file.find_last_of('.'))+ ".out.h5";
        }
    }

}
