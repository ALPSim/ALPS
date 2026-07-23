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

#include "ndim_spin.hpp"

#include <alps/parseargs.hpp>
#include <alps/stop_callback.hpp>
#include <alps/ngs/make_parameters_from_xml.hpp>

#include <boost/chrono.hpp>
#include <boost/filesystem/path.hpp>

#include <string>
#include <iostream>
#include <stdexcept>

typedef ndim_spin_sim<3> sim_type;

int main(int argc, char *argv[]) {
    try {
        // read commandline args
        alps::parseargs options(argc, argv);

        // name the checkpoint file after the input file
        std::string checkpoint_file = options.input_file.substr(0, options.input_file.find_last_of('.')) +  ".checkpt.h5";

        // read in parameters from the input file, can be xml, hd5f, or plain text
        alps::parameters_type<sim_type>::type parameters;
        if (boost::filesystem::path(options.input_file).extension().string() == ".xml") {
            parameters = alps::make_parameters_from_xml(options.input_file);
        }
        else if (boost::filesystem::path(options.input_file).extension().string() == ".h5") {
            alps::hdf5::archive(options.input_file)["/parameters"] >> parameters;
        }
        else {
            parameters = alps::parameters_type<sim_type>::type(options.input_file);
        }

        // initialize the simulation
        sim_type sim(parameters);

        // if the -c option is specified, load state from existing checkpoint file
        if (options.resume)
            sim.load(checkpoint_file);

        // run the simulation with default or commandline given time limit and save the state
        sim.run(alps::stop_callback(options.timelimit));
        sim.save(checkpoint_file);

        // get the results and write an overeview to stdout
        using alps::collect_results;
        alps::results_type<sim_type>::type results = collect_results(sim);
        std::cout << results << std::endl;

        // store parameters and results to the output file (default named after input file)
        alps::hdf5::archive ar(options.output_file, "w");
        ar["/parameters"] << parameters;
        ar["/simulation/results"] << results;

    } catch (std::exception const & e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
