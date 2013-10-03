/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * This software is part of the ALPS libraries, published under the ALPS           *
 * Library License; you can use, redistribute it and/or modify it under            *
 * the terms of the license, either version 1 or (at your option) any later        *
 * version.                                                                        *
 *                                                                                 *
 * You should have received a copy of the ALPS Library License along with          *
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also       *
 * available from http://alps.comp-phys.org/.                                      *
 *                                                                                 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        *
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT       *
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE       *
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,     *
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER     *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "ising.hpp"

#include <alps/ngs.hpp>
#include <alps/ngs/scheduler/mpi_adapter.hpp>
#include <alps/ngs/make_parameters_from_xml.hpp>

#include <boost/chrono.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>

#include <string>
#include <iostream>
#include <stdexcept>

typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

int main(int argc, char *argv[]) {

    try {
        boost::mpi::environment env(argc, argv);
        boost::mpi::communicator comm_world;

        alps::parseargs options(argc, argv);

        std::vector<std::string> input_files;

        tokenizer input_file_list(options.input_file, boost::char_separator<char>(","));
        for (tokenizer::const_iterator it = input_file_list.begin(); it != input_file_list.end(); ++it)
            input_files.push_back(*it);

        std::size_t color = comm_world.rank() * input_files.size() / comm_world.size();
        boost::mpi::communicator comm_local = comm_world.split(color);

        // boost::mpi::communicator & comm_local = comm_world;

        std::string checkpoint_file = input_files[color].substr(0, input_files[color].find_last_of('.')) 
                                    +  ".clone" + boost::lexical_cast<std::string>(comm_local.rank()) + ".h5";

        alps::parameters_type<ising_sim>::type parameters;
        if (comm_local.rank() > 0);
        else if (boost::filesystem::extension(input_files[color]) == ".xml")
            parameters = alps::make_parameters_from_xml(options.input_file);
        else if (boost::filesystem::extension(input_files[color]) == ".h5")
            alps::hdf5::archive(input_files[color])["/parameters"] >> parameters;
        else
            parameters = alps::parameters_type<ising_sim>::type(input_files[color]);
        broadcast(comm_local, parameters);

        alps::mpi_adapter<ising_sim> sim(parameters, comm_local, alps::check_schedule(options.tmin, options.tmax));

        if (options.resume)
            sim.load(checkpoint_file);

        // TODO: how do we handle signels in mpi context? do we want to handle these in the callback or in the simulation?
        // do not use stop_callback_mpi: we do not want an bcast after every sweep!
        sim.run(alps::stop_callback(comm_local, options.timelimit));

        sim.save(checkpoint_file);

        using alps::collect_results;
        alps::results_type<ising_sim>::type results = collect_results(sim);

        if (comm_local.rank() == 0) {
            std::cout << results << std::endl;
            alps::hdf5::archive ar(options.output_file, "w");
            ar["/parameters"] << parameters;
            ar["/simulation/results"] << results;
        }

    } catch (std::exception const & e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
