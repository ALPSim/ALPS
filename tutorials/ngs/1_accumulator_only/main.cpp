/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "ising.hpp"

#include <alps/ngs.hpp>

#include <boost/chrono.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>

#include <string>
#include <iostream>
#include <stdexcept>

struct stop_callback {
    stop_callback(std::size_t timelimit)
        : limit(timelimit)
        , start(boost::chrono::high_resolution_clock::now())
    {}

    bool operator()() {
        return !signals.empty() 
            || (limit.count() > 0 && boost::chrono::high_resolution_clock::now() > start + limit);
    }

    boost::chrono::duration<std::size_t> limit;
    alps::ngs::signal signals;
    boost::chrono::high_resolution_clock::time_point start;
};

int main(int argc, char *argv[]) {

    if (argc != 3) {
        std::cerr << "Usage: timelimit parameter-file" << std::endl;
        return EXIT_FAILURE;
    }

    try {

        std::size_t timelimit = boost::lexical_cast<std::size_t>(argv[1]);
        boost::filesystem::path input_file = argv[2];
        std::string basename = std::string(argv[2]).substr(0, std::string(argv[2]).find_last_of('.'));
        boost::filesystem::path checkpoint_file = basename + ".clone0.h5";
        boost::filesystem::path output_file = basename +  ".out.h5";

        alps::parameters_type<ising_sim>::type parameters(input_file);

        ising_sim sim(parameters);

        sim.run(stop_callback(timelimit));
        // TODO: either drop checkpoint file or add support for resuming (in a second example?)
        sim.save(checkpoint_file);

        using alps::collect_results;
        alps::results_type<ising_sim>::type results = collect_results(sim);

        std::cout << results << std::endl;
        alps::hdf5::archive ar(output_file, "w");
        ar["/parameters"] << parameters;
        ar["/simulation/results"] << results;

    } catch (std::exception const & e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
