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

typedef ising_simulation<alps::mcthreadedsim<alps::mcbase> > Simulation;

int main(int argc, char *argv[]) {

    // load command line argumen
    alps::mcoptions options(argc, argv);

    // load parameterfile
    alps::parameters_type<Simulation>::type params(options.input_file);


    // create simulation
    Simulation s(params);

    // resume if --continue is passed
    if (options.resume)
        s.load(params.value_or_default("DUMP", "dump").str());

    // make callback
    alps::threaded_callback_wrapper stopper(boost::bind<bool>(&alps::basic_stop_callback, options.time_limit));

    // runs simulation
    boost::thread thread(boost::bind<bool>(&Simulation::run, boost::ref(s), stopper));
    // output limits
    boost::posix_time::ptime progress_time = boost::posix_time::second_clock::local_time();
    boost::posix_time::ptime checkpoint_time = boost::posix_time::second_clock::local_time();
    do {
        usleep(0.1 * 1e6);
        // print progress every 5s
        if (boost::posix_time::second_clock::local_time() > progress_time + boost::posix_time::seconds(5)) {
            std::cout << "progress: " << s.fraction_completed() << std::endl;
            progress_time = boost::posix_time::second_clock::local_time();
        }
        // checkpoint every 15 s
        if (boost::posix_time::second_clock::local_time() > checkpoint_time + boost::posix_time::seconds(10)) {
            std::cout << "checkpointing ... " << std::endl;
            // save observables to hdf5
            s.save(params.value_or_default("DUMP", "dump").str());
            checkpoint_time = boost::posix_time::second_clock::local_time();
        }
    } while (!stopper.check());
    thread.join();

    // save observables to hdf5
    s.save(params.value_or_default("DUMP", "dump").str());

    alps::results_type<Simulation>::type results = collect_results(s);
    {
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
    }

    // save results to hdf5
    save_results(results, params, options.output_file, "/simulation/results");

}

