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

// TODO: merge changes back from main_controlthread_tiny!

#include "ising.hpp"

using namespace alps;

typedef ising_sim<multithread<base> > sim_type;

// TODO: synchonize with main_thread_tiny
int main(int argc, char *argv[]) {

    mcoptions options(argc, argv);
    parameters_type<sim_type>::type params(hdf5::archive(options.input_file));
    sim_type sim(params);

    if (options.resume)
        sim.load(params["DUMP"] | "dump");

    threaded_callback_wrapper stopper(boost::bind<bool>(&stop_callback, options.time_limit));
    boost::thread thread(&sim_type::run, boost::ref(sim), static_cast<boost::function<bool()> >(boost::ref(stopper)));

  //TODO: create one sensible exampel: progress every 5 minutes, checkpoint every hour
  // and call this a "tiny" example main_threaded_tiny
  
    boost::posix_time::ptime progress_time = boost::posix_time::second_clock::local_time();
    boost::posix_time::ptime checkpoint_time = boost::posix_time::second_clock::local_time();
    do {
        alps::sleep(0.1 * 1e9);

        //TODO: print progress every 5s
        if (boost::posix_time::second_clock::local_time() > progress_time + boost::posix_time::seconds(5)) {
            std::cout << "progress: " << sim.fraction_completed() << std::endl;
            progress_time = boost::posix_time::second_clock::local_time();
        }

        //TODO: checkpoint every 15 s
        if (boost::posix_time::second_clock::local_time() > checkpoint_time + boost::posix_time::seconds(15)) {
            std::cout << "checkpointing ... " << std::endl;
            sim.save(params["DUMP"] | "dump");
            checkpoint_time = boost::posix_time::second_clock::local_time();
        }

    } while (!stopper.check());
    thread.join();

    sim.save(params["DUMP"] | "dump");

    results_type<sim_type>::type results = collect_results(sim);

    save_results(results, params, options.output_file, "/simulation/results");

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


}
