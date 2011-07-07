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

int main(int argc, char *argv[]) {
/*
    // load command line argumen
    alps::mcoptions options(argc, argv);

    // initialize mpi
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator c;
    
    // load parameterfile
    alps::parameters_type<alps::mcmpisim<alps::mcthreadedsim<ising_simulation> > >::type params(options.input_file);

    // create simulation
    alps::mcmpisim<alps::mcthreadedsim<ising_simulation> > s(params, c);
    
    // resume if --continue is passed
    if (options.resume)
        s.load(params.value_or_default("DUMP", "dump").str() + boost::lexical_cast<std::string>(c.rank()));

    // runs simulation
    s.run(boost::bind(&alps::basic_stop_callback, boost::posix_time::second_clock::local_time(), options.time_limit));

    // save observables to hdf5
    s.save(params.value_or_default("DUMP", "dump").str() + boost::lexical_cast<std::string>(c.rank()));

    alps::results_type<alps::mcmpisim<alps::mcthreadedsim<ising_simulation> > >::type results = collect_results(s);
    if (!c.rank()) {
        using namespace alps;
        // print whole result
        std::cout << "Correlations:           " << short_print(results["Correlations"]) << std::endl;
        std::cout << "Sign:                   " << short_print(results["Sign"]) << std::endl;

        // access meand and error
        std::cout << "Mean of Energy:         " << short_print(results["Energy"].mean<double>()) << std::endl;
        std::cout << "Error of Energy:        " << short_print(results["Energy"].mean<double>()) << std::endl;
        std::cout << "Mean of Correlations:   " << short_print(results["Correlations"].mean<std::vector<double> >()) << std::endl;

        // results can be used with numeric operators
        std::cout << "2 * Energy / 13:        " << short_print(2. * results["Energy"] / 13.) << std::endl;
        std::cout << "1 / Correlations * Sign " << short_print(1. / results["Correlations"] * results["Sign"]) << std::endl;
        std::cout << "-2 * Energy:            " << short_print(-2. * results["Energy"]) << std::endl;
        std::cout << "Correlations / Sign:    " << short_print(results["Correlations"] / results["Sign"]) << std::endl;
        std::cout << "Sign / Correlations:    " << short_print(results["Sign"] / results["Correlations"]) << std::endl;
        
        // more complex functions can be applied
        std::cout << "Sin(Energy):            " << short_print(sin(results["Energy"])) << std::endl;
        std::cout << "Tanh(Correlations):     " << short_print(tanh(results["Correlations"])) << std::endl;

        // save results to hdf5
        save_results(results, params, options.output_file, "/simulation/results");
    }
*/
}

