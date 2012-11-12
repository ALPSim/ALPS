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

#include "src/ising.hpp"

#include <alps/ngs/scheduler/proto/mpisim.hpp>
#include <alps/ngs/scheduler/proto/controlthreadsim.hpp>

#include <boost/lexical_cast.hpp>

using namespace alps;

typedef mpisim_ng<controlthreadsim_ng<ising_sim> > sim_type;

int main(int argc, char *argv[]) {

    mcoptions options(argc, argv);
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator c;

    parameters_type<sim_type>::type params;
    if (c.rank() == 0)
        hdf5::archive(options.input_file)["/parameters"] >> params;
    broadcast(c, params);

    sim_type sim(params, c);

    if (options.resume)
        sim.load(params["DUMP"] | "checkpoint");

    sim.run(boost::bind(&stop_callback, options.time_limit));

    boost::posix_time::ptime progress_time = boost::posix_time::second_clock::local_time();
    boost::posix_time::ptime checkpoint_time = boost::posix_time::second_clock::local_time();
    do {

        alps::sleep(0.1 * 1e9);

        if (c.rank() == 0) {
            if (boost::posix_time::second_clock::local_time() > progress_time + boost::posix_time::seconds(5)) {
                std::cout << "progress: " << std::setprecision(2) << 100 * sim.fraction_completed() << "%" << std::endl;
                progress_time = boost::posix_time::second_clock::local_time();
            }

            if (boost::posix_time::second_clock::local_time() > checkpoint_time + boost::posix_time::seconds(13)) {
                std::cout << "checkpointing rank " << c.rank() << "... " << std::endl;
                sim.save(params["DUMP"] | "checkpoint");
                checkpoint_time = boost::posix_time::second_clock::local_time();
            }
        }

        sim.check_communication();

    } while (sim.status() != ising_sim::finished);

    sim.save((params["DUMP"] | "checkpoint"));

    results_type<sim_type>::type results = collect_results(sim);

    if (c.rank() == 0) {
        std::cout << results << std::endl;
        save_results(results, params, options.output_file, "/simulation/results");
    }
}
