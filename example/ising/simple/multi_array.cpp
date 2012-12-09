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

#include <alps/ngs/api.hpp>
#include <alps/ngs/hdf5.hpp>
#include <alps/ngs/config.hpp>
#include <alps/ngs/params.hpp>
#include <alps/ngs/alea/accumulator_set.hpp>
#include <alps/ngs/make_parameters_from_xml.hpp>

#include <alps/random/mersenne_twister.hpp>

#include <boost/chrono.hpp>
#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>

#ifndef ALPS_NGS_USE_NEW_ALEA
#error "this example only works with the new alea"
#endif

class multi_array_sim {

    public:
        
        typedef alps::params parameters_type;
        typedef std::map<std::string, std::pair<unsigned long, alps::multi_array<double, 2> > > results_type;
        typedef std::vector<std::string> result_names_type;

        multi_array_sim(parameters_type const & parameters)
            : params(parameters)
            , random(boost::mt19937((parameters["SEED"] | 42)), boost::uniform_real<>())
            , steps(0)
        {
            using namespace alps::alea;
            measurements << make_accumulator(
                  "M"
                , detail::accumulator_wrapper(accumulator<alps::multi_array<double, 2>, features<tag::mean> >())
            );
        }

        void update() {};

        void measure() {
            static const unsigned size = 100;
            ++steps;
            alps::multi_array<double, 2> values(size, size);
            std::generate(values.data(), values.data() + size * size, random);
            measurements["M"] << values;
        };
        
        double fraction_completed() const {
            static const double total = 100.;
            return steps / total;
        }

        void save(boost::filesystem::path const & filename) const {
            alps::hdf5::archive ar(filename, "w");
            ar << *this;
        }

        void load(boost::filesystem::path const & filename) {
            alps::hdf5::archive ar(filename);
            ar >> *this;
        }

        bool run(boost::function<bool ()> const & stop_callback) {
            bool stopped = false;
            do {
                update();
                measure();
            } while(!(stopped = stop_callback()) && fraction_completed() < 1.);
            return !stopped;
        }
        
        result_names_type result_names() const {
            result_names_type names;
            for(alps::alea::accumulator_set::const_iterator it = measurements.begin(); it != measurements.end(); ++it)
                names.push_back(it->first);
            return names;
        }

        result_names_type unsaved_result_names() const {
            return result_names_type(); 
        }

        results_type collect_results() const {
            return collect_results(result_names());
        }

        results_type collect_results(result_names_type const & names) const {
            results_type partial_results;
            for(result_names_type::const_iterator it = names.begin(); it != names.end(); ++it) {
                lps::multi_array<double, 2> sum = measurements[*it].get<alps::multi_array<double, 2> >().mean();
                sum *= measurements[*it].count();
                partial_results[*it] = sum;
            }
            return partial_results;
        }

        void save(alps::hdf5::archive & ar) const {
            ar["/parameters"] << params;

            std::string context = ar.get_context();
            ar.set_context("/simulation/realizations/0/clones/0");
            ar["measurements"] << measurements;

            ar.set_context("checkpoint");
            ar["steps"] << steps;

            {
                std::ostringstream os;
                os << random.engine();
                ar["engine"] << os.str();
            }

            ar.set_context(context);
        }

        void load(alps::hdf5::archive & ar) {
            ar["/parameters"] >> params;

            std::string context = ar.get_context();
            ar.set_context("/simulation/realizations/0/clones/0");
            ar["measurements"] >> measurements;

            ar.set_context("checkpoint");
            ar["steps"] >> steps;

            {
                std::string state;
                ar["engine"] >> state;
                std::istringstream is(state);
                is >> random.engine();
            }

            ar.set_context(context);
        }

    private:

        parameters_type params;
        boost::variate_generator<boost::mt19937, boost::uniform_real<> > mutable random;
        alps::alea::accumulator_set measurements;
    
        int steps;
};

bool stop_callback() {
    return true;
}

int main(int argc, char *argv[]) {

    try {
        alps::parameters_type<multi_array_sim>::type parameters;
        parameters["SEED"] = 66;

        multi_array_sim sim(parameters);

        sim.run(&stop_callback);
        
        sim.save("checkpoint.h5");
        
        using alps::collect_results;
        alps::results_type<multi_array_sim>::type results = collect_results(sim);

        alps::hdf5::archive ar("result.h5", "w");
        ar["/parameters"] << parameters;
        ar["/simulation/results"] << results;

    } catch (std::exception const & e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Caught unknown exception" << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
