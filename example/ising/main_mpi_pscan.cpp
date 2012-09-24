/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2012 by Lukas Gamper <gamperl@gmail.com>                   *
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

#include <boost/lexical_cast.hpp>

using namespace alps;

typedef mcmpisim<ising_sim<mcbase> > sim_type;

int main(int argc, char *argv[]) {

    mcoptions options(argc, argv);
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm_world;

    std::size_t color = comm_world.rank() % 2;

    boost::mpi::communicator comm_local = comm_world.split(color);

    parameters_type<sim_type>::type params;
    if (!comm_local.rank()) {
        hdf5::archive ar(options.input_file + "." + boost::lexical_cast<std::string>(color));
        ar["/parameters"] >> params;
    }
    broadcast(comm_local, params);
    
    std::cout << boost::lexical_cast<std::string>(comm_world.rank()) + " " + boost::lexical_cast<std::string>(comm_local.rank()) + " " + params["L"] << std::endl;

    sim_type sim(params, comm_local);

    if (options.resume)
        sim.load(
              (params["DUMP"] | "dump")
            + "." + boost::lexical_cast<std::string>(color) 
            + "." + boost::lexical_cast<std::string>(comm_local.rank())
        );

    sim.run(boost::bind(&basic_stop_callback, options.time_limit));

    sim.save(
          (params["DUMP"] | "dump")
        + "." + boost::lexical_cast<std::string>(color) 
        + "." + boost::lexical_cast<std::string>(comm_local.rank())
    );

    results_type<sim_type>::type results = collect_results(sim);

    if (!comm_local.rank()) {
        std::cout << "Correlations:           " << results["Correlations"] << std::endl;
        std::cout << "Energy:                 " << results["Energy"] << std::endl;

        std::cout << "Mean of Energy:         " << results["Energy"].mean<double>() << std::endl;
        std::cout << "Error of Energy:        " << results["Energy"].error<double>() << std::endl;
        std::cout << "Mean of Correlations:   " << short_print(results["Correlations"].mean<std::vector<double> >()) << std::endl;

        std::cout << "-2 * Energy / 13:       " << -2. * results["Energy"] / 13. << std::endl;
        std::cout << "1 / Correlations        " << 1. / results["Correlations"] << std::endl;
        std::cout << "Energy - Magnetization: " << results["Energy"] - results["Magnetization"] << std::endl;

        std::cout << "Sin(Energy):            " << sin(results["Energy"]) << std::endl;
        std::cout << "Tanh(Correlations):     " << tanh(results["Correlations"]) << std::endl;

        save_results(results, params, options.output_file + "." + boost::lexical_cast<std::string>(color), "/simulation/results");

    }

}
