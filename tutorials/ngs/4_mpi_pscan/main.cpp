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
#include <alps/mcmpiadapter.hpp>
#include <alps/ngs/make_parameters_from_xml.hpp>

#include <boost/chrono.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem/path.hpp>

#include <string>
#include <iostream>
#include <stdexcept>

// TODO: think about a nice task q as a separate more advanced tutorial

// TODO: why not out the contents of main into a template function, templated on the simulation type? then people can trivially reuse these codes
// and we might also put the ones we need for ALPS applications into the library for use in our codes and mention that here

typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

int main(int argc, char *argv[]) {

    try {
        boost::mpi::environment env(argc, argv);
        boost::mpi::communicator comm_world;

        // TODO: accept space sparated input files, not comma
        alps::parseargs options(argc, argv);

        std::vector<std::string> input_files;

        tokenizer input_file_list(options.input_file, boost::char_separator<char>(","));
        for (tokenizer::const_iterator it = input_file_list.begin(); it != input_file_list.end(); ++it)
            input_files.push_back(*it);
        std::size_t tasks_done = 0;

        while (tasks_done < input_files.size()) {
            // TODO: add comments
            std::size_t color = tasks_done;
            if (comm_world.size() < input_files.size())
                color += comm_world.rank();
            else
                color += comm_world.rank() * input_files.size() / comm_world.size();
            tasks_done += comm_world.size();

            std::cout << color << std::endl;

            boost::mpi::communicator comm_local = comm_world.split(color);
            if (color >= input_files.size())
                break;

            std::string infile = input_files[color];
            std::string checkpoint_file = infile.substr(0, infile.find_last_of('.')) 
                                        +  ".clone" + boost::lexical_cast<std::string>(comm_local.rank()) + ".h5";
          
            alps::parameters_type<ising_sim>::type parameters;
            if (comm_local.rank() > 0);
            else if (boost::filesystem::path(infile).extension().string()== ".xml")
                parameters = alps::make_parameters_from_xml(infile);
            else if (boost::filesystem::path(infile).extension().string() == ".h5")
                alps::hdf5::archive(infile)["/parameters"] >> parameters;
            else
                parameters = alps::parameters_type<ising_sim>::type(infile);
            broadcast(comm_local, parameters);

            alps::mcmpiadapter<ising_sim> sim(parameters, comm_local, alps::check_schedule(options.tmin, options.tmax));

            if (options.resume)
                sim.load(checkpoint_file);

            sim.run(alps::stop_callback(comm_local, options.timelimit));

            sim.save(checkpoint_file);

            using alps::collect_results;
            alps::results_type<ising_sim>::type results = collect_results(sim);

            if (comm_local.rank() == 0) {
                std::cout << results << std::endl;
                std::string output_file = infile.substr(0, infile.find_last_of('.')) + ".out.h5";
                alps::hdf5::archive ar(output_file, "w");
                ar["/parameters"] << parameters;
                ar["/simulation/results"] << results;
            }

        }
    } catch (std::exception const & e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
