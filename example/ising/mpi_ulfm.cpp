/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Donjan Rodic <drodic@phys.ethz.ch>                 *
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

#include <alps/ngs/scheduler/proto/mpisim_ulfm.hpp>
#include <alps/ngs/ulfm.hpp>

#include <boost/lexical_cast.hpp>

using namespace alps;

typedef mpisim_ng<ising_sim> sim_type;

int main(int argc, char *argv[]) {

    mcoptions options(argc, argv);
    boost::mpi::environment env(argc, argv);

    boost::mpi::communicator c = ngs::duplicate_comm(MPI_COMM_WORLD);

    parameters_type<sim_type>::type params;
    if (c.rank() == 0)
        hdf5::archive(options.input_file)["/parameters"] >> params;
    ULFM_REPEAT_IF_RANKS_FAIL_BEGIN
        broadcast(c, params);
    ULFM_REPEAT_IF_RANKS_FAIL_END(c)

    sim_type sim(params, c);

    if (options.resume)
        sim.load(params["DUMP"] | "checkpoint");

    ULFM_REPEAT_IF_RANKS_FAIL_BEGIN
        sim.run(boost::bind(&stop_callback, options.time_limit));
    ULFM_REPEAT_IF_RANKS_FAIL_END(c)

    ULFM_REPEAT_IF_RANKS_FAIL_BEGIN
        sim.save(params["DUMP"] | "checkpoint");
    ULFM_REPEAT_IF_RANKS_FAIL_END(c)

    results_type<sim_type>::type results;
    ULFM_REPEAT_IF_RANKS_FAIL_BEGIN
        results = collect_results(sim);
    ULFM_REPEAT_IF_RANKS_FAIL_END(c)

    if (c.rank() == 0) {
        std::cout << results << std::endl;
        save_results(results, params, options.output_file, "/simulation/results");
    }

}
