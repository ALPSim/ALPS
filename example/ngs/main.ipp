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

bool stop_callback(boost::posix_time::ptime const & end_time) {
    return alps::mcsignal() || boost::posix_time::second_clock::local_time() > end_time;
}
int main(int argc, char *argv[]) {
    alps::mcoptions options(argc, argv);
    if (options.valid && options.type == alps::mcoptions::SINGLE) {
        alps::parameters_type<simulation_type>::type params(options.input_file);
        std::string seed = params["SWEEPS"];
        simulation_type s(params);
        s.run(boost::bind(&stop_callback, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(options.time_limit)));
        s.save_collected(options.output_file);
        alps::results_type<alps::mcmpisim<simulation_type> >::type results = collect_results(s);
        for (alps::results_type<alps::mcmpisim<simulation_type> >::type::const_iterator it = results.begin(); it != results.end(); ++it)
            std::cout << std::fixed << std::setprecision(5) << it->first << ": " << it->second->to_string() << std::endl;
    } else if(options.valid && options.type == alps::mcoptions::MPI) {
        boost::mpi::environment env(argc, argv);
        boost::mpi::communicator c;
        alps::parameters_type<simulation_type>::type params(options.input_file);
        alps::mcmpisim<simulation_type> s(params, c);
        s.run(boost::bind(&stop_callback, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(options.time_limit)));
        s.save_collected(options.output_file);
        alps::results_type<alps::mcmpisim<simulation_type> >::type results = collect_results(s);
        for (alps::results_type<alps::mcmpisim<simulation_type> >::type::const_iterator it = results.begin(); it != results.end(); ++it)
            std::cout << std::fixed << std::setprecision(5) << it->first << ": " << it->second->to_string() << std::endl;
    }
}
