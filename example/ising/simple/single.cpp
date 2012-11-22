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
#include <alps/ngs/hdf5.hpp>
#include <alps/ngs/config.hpp>
#include <alps/ngs/signal.hpp>
#include <alps/ngs/params.hpp>
#include <alps/ngs/mcresults.hpp> // TODO: replace by new alea
#include <alps/ngs/mcobservables.hpp> // TODO: replace by new alea
#include <alps/ngs/make_parameters_from_xml.hpp>

#include <alps/random/mersenne_twister.hpp>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <vector>
#include <string>

class ising_sim {

    public:

        typedef alps::params parameters_type;
        typedef alps::mcresults results_type;
        typedef std::vector<std::string> result_names_type;

        // TODO: use new alea
        #ifdef ALPS_NGS_USE_NEW_ALEA
            typedef alea::accumulator_set observables_type;
        #else
            typedef alps::mcobservables observables_type;
        #endif

        ising_sim(boost::filesystem::path const & filename)
            // TODO: this ist not the best solution - any idea?
            : parameters(alps::make_parameters_from_xml(filename))
            , random(boost::mt19937((parameters["SEED"] | 42)), boost::uniform_real<>())
            , length(parameters["L"])
            , sweeps(0)
            , thermalization_sweeps(int(parameters["THERMALIZATION"]))
            , total_sweeps(int(parameters["SWEEPS"]))
            , beta(1. / double(parameters["T"]))
            , spins(length)
        {
            for(int i = 0; i < length; ++i)
                spins[i] = (random() < 0.5 ? 1 : -1);
            measurements
                << alps::ngs::RealObservable("Energy")
                << alps::ngs::RealObservable("Magnetization")
                << alps::ngs::RealObservable("Magnetization^2")
                << alps::ngs::RealObservable("Magnetization^4")
                << alps::ngs::RealVectorObservable("Correlations")
            ;
            alps::ngs::signal::listen();
        }

        void update() {
            for (int j = 0; j < length; ++j) {
                using std::exp;
                int i = int(double(length) * random());
                int right = ( i + 1 < length ? i + 1 : 0 );
                int left = ( i - 1 < 0 ? length - 1 : i - 1 );
                double p = exp( 2. * beta * spins[i] * ( spins[right] + spins[left] ));
                if ( p >= 1. || random() < p )
                    spins[i] = -spins[i];
            }
        };

        void measure() {
            sweeps++;
            if (sweeps > thermalization_sweeps) {
                double tmag = 0;
                double ten = 0;
                double sign = 1;
                std::vector<double> corr(length);
                for (int i = 0; i < length; ++i) {
                    tmag += spins[i];
                    sign *= spins[i];
                    ten += -spins[i] * spins[ i + 1 < length ? i + 1 : 0 ];
                    for (int d = 0; d < length; ++d)
                        corr[d] += spins[i] * spins[( i + d ) % length ];
                }
                std::transform(corr.begin(), corr.end(), corr.begin(), boost::lambda::_1 / double(length));
                ten /= length;
                tmag /= length;
                measurements["Energy"] << ten;
                measurements["Magnetization"] << tmag;
                measurements["Magnetization^2"] << tmag * tmag;
                measurements["Magnetization^4"] << tmag * tmag * tmag * tmag;
                measurements["Correlations"] << corr;
            }
        };

        double fraction_completed() const {
            return (sweeps < thermalization_sweeps ? 0. : ( sweeps - thermalization_sweeps ) / double(total_sweeps));
        }

        void save(boost::filesystem::path const & filename) const {
            alps::hdf5::archive ar(filename, "w");
            ar["/"] << *this;
        }

        void load(boost::filesystem::path const & filename) {
            alps::hdf5::archive ar(filename);
            ar["/"] >> *this;
        }

       // make these private but hdf5::arhive is a friend
       void save(alps::hdf5::archive & ar) const {
            ar["/parameters"] << parameters;
            ar["/simulation/realizations/0/clones/0/results"] << measurements;
            // TODO: should we save sweeps? or state? yes! of course!
            // TODO: realization and clone number
        }

        // TODO: do we want to load the parameters?
        void load(alps::hdf5::archive & ar) {
            ar["/simulation/realizations/0/clones/0/results"] >> measurements;
        }

        bool run(boost::function<bool ()> const & stop_callback) {
          bool stopped=false;
          do {
                update();
                measure();
            } while(!(stopped=stop_callback()) && fraction_completed() < 1.);
            return !stopped;
        }

        // implement a nice keys(m) function
        result_names_type result_names() const {
            result_names_type names;
            for(observables_type::const_iterator it = measurements.begin(); it != measurements.end(); ++it)
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
            for(result_names_type::const_iterator it = names.begin(); it != names.end(); ++it)
                // TODO: this is ugly make measurements[*it]
                partial_results.insert(*it, alps::mcresult(measurements[*it]));
            return partial_results;
        }

        // TODO: how do we want to call that?
  //    double get_random() const { return random(); }

  // rename to parameters()
        parameters_type const & get_parameters() const { return parameters; }

  //observables_type const & get_measurements() const { return measurements; }

    private:

        parameters_type parameters;
        boost::variate_generator<boost::mt19937, boost::uniform_real<> > mutable random;
        observables_type measurements;

        int length;
        int sweeps;
        int thermalization_sweeps;
        int total_sweeps;
        double beta;
        std::vector<int> spins;
};

int main(int argc, char *argv[]) {

    try {

        // TODO: improve this! how should we specify the options?
      // either explicitly parse the options or use Boost.ProgramOptions
      
      // usage: [-T timelimit] [-c] input_file [output_file]
        alps::mcoptions options(argc, argv);

      alps::parameters parms;
      // load parms from input file
      
      
        ising_sim sim(parms);

      // use '-c checkpointfile' instead of 'inputfile' for continuing from checkpoint
        if (!options.checkpoint_file.empty() && boost::filesystem::exists(options.checkpoint_file))
            sim.load(options.checkpoint_file);

        {
            using alps::stop_callback;
            sim.run(boost::bind(&alps::stop_callback, options.time_limit));
        }

        if (!options.checkpoint_file.empty())
            sim.save(options.checkpoint_file);

        {
            using alps::save_results;
            using alps::collect_results;
            alps::results_type<ising_sim>::type results = collect_results(sim);
            std::cout << results << std::endl;
          // explicitly use an archive: put the function here
            save_results(results, sim.get_parameters(), options.output_file, "/simulation/results");
        }

    } catch (std::exception const & e) {
        std::cerr << "Caught exception: " << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Caught unknown exception" << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
