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

#include "ising.hpp"

using namespace alps;

typedef ising_sim<base> sim_type;

int main(int argc, char *argv[]) {

    mcoptions options(argc, argv);
    parameters_type<sim_type>::type params(hdf5::archive(options.input_file));
    sim_type sim(params);

    if (options.resume)
        sim.load(params.value_or_default("DUMP", "dump"));

    sim.run(boost::bind(&basic_stop_callback, options.time_limit));

    // save simulation to checkpoint
    sim.save(params.value_or_default("DUMP", "dump"));

    results_type<sim_type>::type results = collect_results(sim);

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

    save_results(results, params, options.output_file, "/simulation/results");

}
