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

#include "ising.hpp"

#include <alps/mc/api.hpp>
#include <alps/utilities/signal.hpp>
#include <alps/ngs/make_parameters_from_file.hpp>

#include <boost/lexical_cast.hpp>
#include <filesystem>

#include <chrono>
#include <string>
#include <iostream>
#include <stdexcept>

struct stop_callback {
    stop_callback(std::size_t timelimit)
        : limit(timelimit)
        , start(std::chrono::steady_clock::now())
    {}

    bool operator()() {
        return !signals.empty() 
            || (limit.count() > 0 && std::chrono::steady_clock::now() > start + limit);
    }

    std::chrono::seconds limit;
    alps::signal signals;
    std::chrono::steady_clock::time_point start;
};

int main(int argc, char *argv[]) {

    if (argc != 3) {
        std::cerr << "Usage: timelimit parameter-file" << std::endl;
        return EXIT_FAILURE;
    }

    try {

        std::size_t timelimit = boost::lexical_cast<std::size_t>(argv[1]);
        std::filesystem::path input_file = argv[2];
        std::string basename = std::string(argv[2]).substr(0, std::string(argv[2]).find_last_of('.'));
        std::filesystem::path checkpoint_file = basename + ".clone0.h5";
        std::filesystem::path output_file = basename +  ".out.h5";

        alps::parameters_type<ising_sim>::type parameters = alps::make_parameters_from_file(input_file);
        ising_sim::define_parameters(parameters);
        if (parameters.has_missing(std::cerr))
            return EXIT_FAILURE;

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
