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
#include <alps/ngs/boost_mpi.hpp>
#include <alps/ngs/mcresults.hpp> // TODO: replace by new alea
#include <alps/ngs/mcobservables.hpp> // TODO: replace by new alea
#include <alps/ngs/make_parameters_from_xml.hpp>
#include <alps/ngs/scheduler/check_schedule.hpp>

#include <alps/random/mersenne_twister.hpp>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <vector>
#include <string>

class ising_sim {

    typedef enum { NOOP_TAG, CHECKPOINT_TAG, FRACTION_TAG, STOP_TAG } tag_type;

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

        ising_sim(parameters_type const & p, boost::mpi::communicator const & c, double t_min = 1, double t_max = 600)
            : communicator(c)
            , parameters(p)
            // TODO: this ist not the best solution - any idea?
            , random(boost::mt19937((parameters["SEED"] | 42) + c.rank()), boost::uniform_real<>())
            , schedule(t_min, t_max)
            , binnumber(p["binnumber"] | std::min(128, 2 * c.size()))
            , fraction(0.)
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
            MPI_Errhandler_set(communicator, MPI_ERRORS_RETURN);
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
            return fraction;
        }

        void save(boost::filesystem::path const & filename) const {
            alps::hdf5::archive ar(filename, "w");
            ar["/checkpoint"] << *this;
        }

        void load(boost::filesystem::path const & filename) {
            alps::hdf5::archive ar(filename);
            ar["/checkpoint"] >> *this;
        }

        void save(alps::hdf5::archive & ar) const {
            ar["/parameters"] << parameters;
            ar["/simulation/realizations/0/clones/0/results"] << measurements;
            // TODO: should we save sweeps? or state?
        }

        // TODO: do we want to load the parameters?
        void load(alps::hdf5::archive & ar) {
            ar["/simulation/realizations/0/clones/0/results"] >> measurements;
        }

        bool run(
              boost::function<bool ()> const & stop_callback
        ) {
            bool done = false, stop;
            do {
                update();
                measure();
                if (schedule.pending() || (stop = stop_callback()))
                    done = communicate(stop);
            } while(!done);
            return !stop;
        }

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
            for(result_names_type::const_iterator it = names.begin(); it != names.end(); ++it) {
                alps::mcresult result(measurements[*it]);
                if (result.count())
                    partial_results.insert(*it, result.reduce(communicator, binnumber));
                else
                    partial_results.insert(*it, result);
            }
            return partial_results;
        }

        // TODO: how do we want to call that?
        double get_random() const { return random(); }

        parameters_type & get_parameters() { return parameters; }
        parameters_type const & get_parameters() const { return parameters; }

        observables_type & get_measurements() { return measurements; }
        observables_type const & get_measurements() const { return measurements; }

    private:

        bool communicate(bool stop) {
            double local_fraction = stop
                ? 1.
                : (sweeps < thermalization_sweeps
                    ? 0.
                    : (sweeps - thermalization_sweeps) / double(total_sweeps)
                  )
            ;
            if (!communicator.rank()) {
                boost::mpi::reduce(communicator, local_fraction, fraction, std::plus<double>(), 0);
                boost::mpi::broadcast(communicator, fraction, 0);
            } else {
                reduce(communicator, local_fraction, std::plus<double>(), 0);
                boost::mpi::broadcast(communicator, fraction, 0);
            }
            schedule.update(fraction);
            return fraction >= 1.;
        }

        boost::mpi::communicator communicator;

        parameters_type parameters;
        boost::variate_generator<boost::mt19937, boost::uniform_real<> > mutable random;
        observables_type measurements;
        alps::check_schedule schedule;

        std::size_t binnumber;
        double fraction;

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
        alps::mcoptions options(argc, argv);

        boost::mpi::environment env(argc, argv);
        boost::mpi::communicator c;

        alps::parameters_type<ising_sim>::type parameters;
        if (c.rank() == 0)
            parameters = alps::make_parameters_from_xml(options.input_file);
        broadcast(c, parameters);

        ising_sim sim(parameters, c);

        // TODO: who should add the suffix to the checkpoint file?
        std::string checkpoint_file = options.checkpoint_file + "." + boost::lexical_cast<std::string>(c.rank()) + ".h5";
        if (!options.checkpoint_file.empty() && boost::filesystem::exists(checkpoint_file))
            sim.load(checkpoint_file);

        {
            using alps::stop_callback;
            sim.run(boost::bind(&alps::stop_callback, options.time_limit));
        }

        if (!options.checkpoint_file.empty())
            sim.save(checkpoint_file);

        {
            using alps::save_results;
            using alps::collect_results;
            alps::results_type<ising_sim>::type results = collect_results(sim);

            if (c.rank() == 0) {
                std::cout << results << std::endl;
                save_results(results, parameters, options.output_file, "/simulation/results");
            }
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
