/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Lukas Gamper <gamperl -at- gmail.com>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#include "ising_parallel.hpp"

using namespace alps;

typedef parallel2<parallel_ising_sim<base> > sim_type; // rename parallel to mpi_parallel or mpi_sim

int main(int argc, char *argv[]) {

    mcoptions options(argc, argv);
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator c;

    parameters_type<sim_type>::type params(c);
    if (c.rank() == 0) {
        hdf5::archive ar(options.input_file);
        ar >> make_pvp("/parameters", params);
    }
    params.broadcast(); // HELP!!!!! broadcast(c,params) or similar - since when do we require params to have a broadcast member function?

    sim_type sim(params, c);

    if (options.resume)
        sim.load(params.value_or_default("DUMP", "dump") + "." + boost::lexical_cast<std::string>(c.rank()));

    sim.run(boost::bind(&basic_stop_callback, options.time_limit));

    sim.save(params.value_or_default("DUMP", "dump") + "." + boost::lexical_cast<std::string>(c.rank()));

    results_type<sim_type>::type results = collect_results(sim);

    if (c.rank() == 0) {
        std::cout << "Energy:                 " << results["Energy"] << std::endl;

        std::cout << "Mean of Energy:         " << results["Energy"].mean<double>() << std::endl;
        std::cout << "Error of Energy:        " << results["Energy"].error<double>() << std::endl;
        std::cout << "Energy^2:               " << results["Energy^2"] << std::endl;
        std::cout << "Magnetization           " << results["Magnetization"] << std::endl;
        std::cout << "Binder Cumulant         " << results["Magnetization^2"] * results["Magnetization^2"] / results["Magnetization^4"]<< std::endl;

        save_results(results, params, options.output_file, "/simulation/results");

    }

}