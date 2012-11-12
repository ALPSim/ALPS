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

// TODO: move to ngs.hpp
#include <alps/ngs/scheduler/tcpserver.hpp>
#include <alps/ngs/scheduler/proto/controlthreadsim.hpp>

#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>

using namespace alps;

typedef controlthreadsim_ng<ising_sim> sim_type;

std::string do_checkpoint(sim_type & sim, std::string const & path) {
    sim.save(path);
    return "done";
}

int main(int argc, char *argv[]) {

    mcoptions options(argc, argv);
    parameters_type<sim_type>::type params(hdf5::archive(options.input_file));
    sim_type sim(params);

    if (options.resume)
        sim.load(params["DUMP"] | "checkpoint");

    sim.run(boost::bind(&stop_callback, options.time_limit));

    tcpserver server(params["PORT"] | 2485); // TODO: which port should we take?
    server.add_action("progress", boost::bind(&boost::lexical_cast<std::string, double>, boost::bind(&sim_type::fraction_completed, boost::ref(sim))));
    server.add_action("checkpoint", boost::bind(&do_checkpoint, boost::ref(sim), params["DUMP"] | "checkpoint"));

    boost::posix_time::ptime progress_time = boost::posix_time::second_clock::local_time();
    boost::posix_time::ptime checkpoint_time = boost::posix_time::second_clock::local_time();
    do {
        alps::sleep(0.1 * 1e9);

        if (boost::posix_time::second_clock::local_time() > progress_time + boost::posix_time::seconds(5)) {
            std::cout << "progress: " << std::setprecision(2) << 100 * sim.fraction_completed() << "%" << std::endl;
            progress_time = boost::posix_time::second_clock::local_time();
        }

        server.poll();
    } while (sim.status() != ising_sim::finished);
    server.stop();

    sim.save(params["DUMP"] | "checkpoint");

    results_type<sim_type>::type results = collect_results(sim);
    std::cout << results << std::endl;
    save_results(results, params, options.output_file, "/simulation/results");
}
