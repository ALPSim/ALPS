/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/ngs/mcoptions.hpp>
#include <alps/ngs/stacktrace.hpp>

#include <boost/program_options.hpp>

#include <iostream>
#include <stdexcept>

namespace alps {

    mcoptions::mcoptions(int argc, char* argv[]) : valid(false), resume(false), type(SINGLE) {
        boost::program_options::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("single", "run single process")
            #ifndef ALPS_NGS_SINGLE_THREAD
                ("threaded", "run in multithread environment")
            #endif
            #ifdef ALPS_HAVE_MPI
                ("mpi", "run in parallel using MPI")
            #endif
            ("continue,c", "load simulation from checkpoint")
            ("time-limit,T", boost::program_options::value<std::size_t>(&time_limit)->default_value(0), "time limit for the simulation")
            ("input-file", boost::program_options::value<std::string>(&input_file), "input file in hdf5 format")
            ("output-file", boost::program_options::value<std::string>(&output_file)->default_value("<unspecified>"), "output file in hdf5 format")
            ("checkpoint-file", boost::program_options::value<std::string>(&checkpoint_file)->default_value(""), "checkpoint file in hdf5 format");
        boost::program_options::positional_options_description p;
        p.add("input-file", 1);
        p.add("output-file", 1);
        p.add("checkpoint-file", 1);
        boost::program_options::variables_map vm;
        boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        boost::program_options::notify(vm);
        if (!(valid = !vm.count("help")))
            std::cout << desc << std::endl;
        else if (input_file.empty())
            throw std::invalid_argument("No job file specified" + ALPS_STACKTRACE);
        if (vm.count("threaded") && vm.count("mpi"))
            type = HYBRID;
        else if (vm.count("threaded"))
            #ifdef ALPS_NGS_SINGLE_THREAD
                throw std::logic_error("Not build with multithread support" + ALPS_STACKTRACE);
            #else
                type = THREADED;
            #endif
        else if (vm.count("mpi")) {
            type = MPI;
            #ifndef ALPS_HAVE_MPI
                throw std::logic_error("Not build with MPI" + ALPS_STACKTRACE);
            #endif
        }
        if (vm.count("continue"))
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
