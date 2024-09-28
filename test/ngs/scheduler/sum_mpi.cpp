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

#include <alps/ngs.hpp>
#include <alps/mcbase.hpp>
#include <alps/parseargs.hpp>
#include <alps/mcmpiadapter.hpp>
#include <alps/stop_callback.hpp>
#include <alps/ngs/make_parameters_from_xml.hpp>

#include <boost/lambda/lambda.hpp>

// Simulation to measure e^(-x*x)
class my_sim_type : public alps::mcbase {

    public:

        my_sim_type(parameters_type const & params, std::size_t seed_offset = 42)
            : alps::mcbase(params, seed_offset)
            , total_count(params["COUNT"])

        {
            measurements << alps::accumulator::RealObservable("SValue")
                         << alps::accumulator::RealVectorObservable("VValue");
        }

        // if not compiled with mpi boost::mpi::communicator does not exists, 
        // so template the function
        template <typename Arg> my_sim_type(parameters_type const & params, Arg comm)
            : alps::mcbase(params, comm)
            , total_count(params["COUNT"])
        {
            measurements << alps::accumulator::RealObservable("SValue")
                         << alps::accumulator::RealVectorObservable("VValue");
        }

        // do the calculation in this function
        void update() {
            double x = random();
            value = exp(-x * x);
        };

        // do the measurements here
        void measure() {
            ++count;
            measurements["SValue"] << value;
            measurements["VValue"] << std::vector<double>(3, value);
        };

        double fraction_completed() const {
            return count / double(total_count);
        }

    private:
        int count;
        int total_count;
        double value;
};

int main(int argc, char *argv[]) {

    try {

        alps::parseargs options(argc, argv);
        boost::mpi::environment env(argc, argv);
        boost::mpi::communicator c;

        alps::parameters_type<my_sim_type>::type params;
        if (c.rank() > 0)
          /* do nothing*/ ;
        else if (boost::filesystem::path(options.input_file).extension().string() == ".xml")
            params = alps::make_parameters_from_xml(options.input_file);
        else if (boost::filesystem::path(options.input_file).extension().string() == ".h5")
            alps::hdf5::archive(options.input_file)["/parameters"] >> params;
        else
            params = alps::parameters_type<my_sim_type>::type(options.input_file);
        broadcast(c, params);

        alps::mcmpiadapter<my_sim_type> my_sim(params, c, alps::check_schedule(options.tmin, options.tmax)); // creat a simulation

        my_sim.run(alps::stop_callback(c, options.timelimit)); // run the simulation

        using alps::collect_results;

        if (c.rank() == 0) { // print the results and save it to hdf5
            alps::results_type<alps::mcmpiadapter<my_sim_type> >::type results = collect_results(my_sim);
            std::cout << "e^(-x*x): " << results["SValue"] << std::endl;
            std::cout << "e^(-x*x): " << results["VValue"] << std::endl;
            using std::sin;
            std::cout << results["SValue"] + 1 << std::endl;
            std::cout << results["SValue"] + results["SValue"] << std::endl;
            // std::cout << results["SValue"] << " " << results["SValue"] * results["SValue"] << sin(results["SValue"]) << std::endl;
            // std::cout << results["SValue"] * results["SValue"] << std::endl;
            // std::cout << << 2. * results["SValue"] / 2. << std::endl;
            // std::cout << sin(results["SValue"]) << std::endl;
            // std::cout << results["VValue"] << " " << 2. * results["VValue"] / 2. << sin(results["SValue"]) << std::endl;
            save_results(results, params, options.output_file, "/simulation/results");
        } else
            collect_results(my_sim);

    } catch(std::exception & ex) {
        std::cerr << ex.what() << std::endl;
        return -1;
    } catch(...) {
        std::cerr << "Fatal Error: Unknown Exception!\n";
        return -2;
    }
}
