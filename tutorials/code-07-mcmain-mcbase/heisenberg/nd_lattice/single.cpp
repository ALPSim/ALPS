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

#include "heisenberg.hpp"

#include <alps/parseargs.hpp>
#include <alps/stop_callback.hpp>
#include <alps/ngs/make_parameters_from_xml.hpp>

#include <boost/chrono.hpp>
#include <boost/filesystem/path.hpp>

#include <string>
#include <iostream>
#include <stdexcept>

int main(int argc, char *argv[]) {

    try {
        alps::parseargs options(argc, argv);
        std::string checkpoint_file = options.input_file.substr(0, options.input_file.find_last_of('.')) +  ".clone0.h5";

        // TODO: make load_params
        alps::parameters_type<heisenberg_sim>::type parameters;
        // TODO: better check the first few bytes. provide an ALPS function to do so
        if (boost::filesystem::path(options.input_file).extension().string() == ".xml") {
            parameters = alps::make_parameters_from_xml(options.input_file);
        }
        else if (boost::filesystem::path(options.input_file).extension().string() == ".h5") {
            alps::hdf5::archive(options.input_file)["/parameters"] >> parameters;
        }
        else {
            parameters = alps::parameters_type<heisenberg_sim>::type(options.input_file);
        }
        //~ alps::Parameters old_parameters = alps::make_deprecated_parameters(parameters);

        //~ heisenberg_sim sim(parameters, old_parameters);
        heisenberg_sim sim(parameters);

        if (options.resume)
            sim.load(checkpoint_file);

        sim.run(alps::stop_callback(options.timelimit));

        sim.save(checkpoint_file);

        using alps::collect_results;
        alps::results_type<heisenberg_sim>::type results = collect_results(sim);

        std::cout << results << std::endl;
        alps::hdf5::archive ar(options.output_file, "w");
        ar["/parameters"] << parameters;
        ar["/simulation/results"] << results;

    } catch (std::exception const & e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
