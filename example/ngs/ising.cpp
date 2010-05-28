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

# include "ising.hpp"

//#define SINGLE_RUN
//#define THREAD_RUN
#define MPI_RUN

#ifdef SINGLE_RUN
    typedef ising_sim sim_type;
#endif
#ifdef THREAD_RUN
    typedef alps::mcthreadsim<ising_sim> sim_type;
#endif
#ifdef MPI_RUN
    typedef alps::mcmpisim<ising_sim> sim_type;
#endif
bool stop_callback(boost::posix_time::ptime const & end_time) {
    return alps::mcsignal() || boost::posix_time::second_clock::local_time() > end_time;
}
int main(int argc, char *argv[]) {
    alps::mcoptions options(argc, argv);
    if (options.valid) {
        sim_type::parameters_type params(options.input_file);
#ifdef MPI_RUN
        boost::mpi::environment env(argc, argv);
        boost::mpi::communicator c;
        sim_type s(params, c);
#else
        sim_type s(params);
#endif
        s.run(boost::bind(&stop_callback, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(options.time_limit)));
#ifdef MPI_RUN
        s.save("sim-" + boost::lexical_cast<std::string>(c.rank()));
        alps::results_type<sim_type>::type results = s.collect_local_results();
        for (alps::results_type<sim_type>::type::const_iterator it = results.begin(); it != results.end(); ++it)
            std::cout << std::fixed << std::setprecision(5) << it->first << " (" << c.rank() << "): " << it->second->to_string() << std::endl;
        if (c.rank()==0) {
            {
                alps::results_type<sim_type>::type results = collect_results(s);
                alps::hdf5::oarchive ar(options.output_file);
                ar << alps::make_pvp("/parameters", params);
                for (alps::results_type<sim_type>::type::const_iterator it = results.begin(); it != results.end(); ++it)
                    ar << alps::make_pvp("/simulation/results/" + it->first, *(it->second));
                for (alps::results_type<sim_type>::type::const_iterator it = results.begin(); it != results.end(); ++it)
                    std::cout << std::fixed << std::setprecision(5) << it->first << ": " << it->second->to_string() << std::endl;
            }
            s.terminate();
        } else
            s.process_requests();
#else
        s.save("sim");
        collect_results(s);
        std::cout << collect_results(s, "Magnetization");
#endif
    }
}
