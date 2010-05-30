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
void run_single(alps::mcoptions const & options) {
    clone_type::parameters_type params(options.input_file);
    clone_type s(params);
    s.run(boost::bind(&stop_callback, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(options.time_limit)));
//    s.save(options.output_file);
    collect_results(s);
    std::cout << collect_results(s, "Magnetization");
}
void run_threaded(alps::mcoptions const & options) {
    alps::mcthreadsim<clone_type>::parameters_type params(options.input_file);
    alps::mcthreadsim<clone_type> s(params);
    s.run(boost::bind(&stop_callback, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(options.time_limit)));
//    s.save(options.output_file);
    collect_results(s);
    std::cout << collect_results(s, "Magnetization");
}
void run_mpi(alps::mcoptions const & options, int argc, char *argv[]) {
    alps::mcmpisim<clone_type>::parameters_type params(options.input_file);
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator c;
    alps::mcmpisim<clone_type> s(params, c);

    s.run(boost::bind(&stop_callback, boost::posix_time::second_clock::local_time() + boost::posix_time::seconds(options.time_limit)));
//    s.save("sim-" + boost::lexical_cast<std::string>(c.rank()));
    alps::results_type<clone_type>::type results = s.collect_local_results();
    for (alps::results_type<clone_type>::type::const_iterator it = results.begin(); it != results.end(); ++it)
        std::cout << std::fixed << std::setprecision(5) << it->first << " (" << c.rank() << "): " << it->second->to_string() << std::endl;
    if (c.rank()==0) {
        {
            alps::results_type<clone_type>::type results = collect_results(s);
            alps::hdf5::oarchive ar(options.output_file);
            ar << alps::make_pvp("/parameters", params);
            for (alps::results_type<clone_type>::type::const_iterator it = results.begin(); it != results.end(); ++it)
                ar << alps::make_pvp("/simulation/results/" + it->first, *(it->second));
            std::cout << results;
        }
        s.terminate();
    } else
        s.process_requests();
}
int main(int argc, char *argv[]) {
    alps::mcoptions options(argc, argv);
    if (options.valid)
        switch(options.type) {
            case alps::mcoptions::SINGLE: run_single(options); break;
            case alps::mcoptions::THREADED: run_threaded(options); break;
            case alps::mcoptions::MPI: run_mpi(options, argc, argv); break;
        }
}