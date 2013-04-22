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
#include <alps/ngs/signal.hpp>
#include <alps/ngs/params.hpp>
#include <alps/ngs/boost_mpi.hpp>
#include <alps/ngs/mcresults.hpp>
#include <alps/ngs/mcobservables.hpp>

#ifdef ALPS_NGS_USE_NEW_ALEA
    #include <alps/ngs/alea/accumulator_set.hpp>
#endif
 
#include <alps/ngs/observablewrappers.hpp>
#include <alps/ngs/make_parameters_from_xml.hpp>
#include <alps/ngs/scheduler/check_schedule.hpp>

#include <alps/random/mersenne_twister.hpp> // TODO: why do we need this?

#include <boost/chrono.hpp>
#include <boost/function.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>

class ising_sim {

    #ifdef ALPS_NGS_USE_NEW_ALEA
        typedef alps::accumulator::accumulator_set observables_type;
    #else
        typedef alps::mcobservables observables_type;
    #endif

    public:

        typedef alps::params parameters_type;
        typedef alps::mcresults results_type;
        typedef std::vector<std::string> result_names_type;

        ising_sim(parameters_type const & parameters, boost::mpi::communicator const & comm, double t_min = 1, double t_max = 600)
            : communicator(comm)
            , params(parameters)
            // TODO: this ist not the nicest solution - any idea?
            , random(boost::mt19937((parameters["SEED"] | 42) + comm.rank()), boost::uniform_real<>())
            , schedule(t_min, t_max)
            , clone(comm.rank())
            , binnumber(parameters["BINNUMBER"] | std::min(128, 2 * comm.size()))
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
        
        //jan: unused, remove?
        double fraction_completed() const {
            return fraction;
        }

        void save(boost::filesystem::path const & filename) const {
            alps::hdf5::archive ar(filename, "w");
            ar["/"] << *this;
        }

        void load(boost::filesystem::path const & filename) {
            alps::hdf5::archive ar(filename);
            ar["/"] >> *this;
        }

        bool run(boost::function<bool ()> const & stop_callback) {
            bool done = false, stopped = false;
            do {
                update();
                measure();
                if ((stopped = stop_callback()) || schedule.pending()) {
                    double local_fraction = (stopped
                        ? 1.
                        : (sweeps < thermalization_sweeps
                            ? 0.
                            : (sweeps - thermalization_sweeps) / double(total_sweeps)
                          )
                    );
                    schedule.update(fraction = boost::mpi::all_reduce(communicator, local_fraction, std::plus<double>()));
                    done = fraction >= 1.;
                }
            } while(!done);
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
            for(result_names_type::const_iterator it = names.begin(); it != names.end(); ++it) {
                alps::mcresult result(measurements[*it]);
                if (result.count())
                    partial_results.insert(*it, result.reduce(communicator, binnumber));
                else
                    partial_results.insert(*it, result);
            }
            return partial_results;
        }

        parameters_type const & parameters() const { return params; }

        void save(alps::hdf5::archive & ar) const {
            ar["/parameters"] << params;
            std::string context = ar.get_context();
            ar.set_context("/simulation/realizations/0/clones/" + boost::lexical_cast<std::string>(clone));
            ar["measurements"] << measurements;

            ar.set_context("checkpoint");
            ar["length"] << length;
            ar["sweeps"] << sweeps;
            ar["thermalization_sweeps"] << thermalization_sweeps;
            ar["beta"] << beta;
            ar["spins"] << spins;
            
            {
                std::ostringstream os;
                os << random.engine();
                ar["engine"] << os.str();
            }

            ar.set_context(context);
        }

        void load(alps::hdf5::archive & ar) {
            ar["/parameters"] >> params; // TODO: do we want to load the parameters?

            std::string context = ar.get_context();
            ar.set_context("/simulation/realizations/0/clones/" + boost::lexical_cast<std::string>(clone));
            ar["measurements"] >> measurements;

            ar.set_context("checkpoint");
            ar["length"] >> length;
            ar["sweeps"] >> sweeps;
            ar["thermalization_sweeps"] >> thermalization_sweeps;
            ar["beta"] >> beta;
            ar["spins"] >> spins;

            {
                std::string state;
                ar["engine"] >> state;
                std::istringstream is(state);
                is >> random.engine();
            }

            ar.set_context(context);
        }

    private:

        boost::mpi::communicator communicator;

        parameters_type params;
        boost::variate_generator<boost::mt19937, boost::uniform_real<> > mutable random;
        observables_type measurements;

        alps::check_schedule schedule;
        double fraction;

        int clone;

        std::size_t binnumber;

        int length;
        int sweeps;
        int thermalization_sweeps;
        int total_sweeps;
        double beta;
        std::vector<int> spins;
};

struct args {
    args(int argc, char *argv[], int rank) {
        boost::program_options::options_description options("Options");
        options.add_options()
            ("continue,c", "load simulation from checkpoint")
            ("timelimit,T", boost::program_options::value<std::size_t>(&timelimit)->default_value(0), "time limit for the simulation")
            ("inputfile", boost::program_options::value<std::string>(&inputfile), "input file in hdf5 or xml format")
            ("outputfile", boost::program_options::value<std::string>(&outputfile)->default_value(""), "output file in hdf5 format")
        ;
        boost::program_options::positional_options_description positional;
        positional
            .add("inputfile", 1)
            .add("outputfile", 1)
        ;

        try {
            boost::program_options::variables_map variables;
            boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(options).positional(positional).run(), variables);
            boost::program_options::notify(variables);

            resume = variables.count("continue");
            checkpointfile = (outputfile.empty()
                ? inputfile.substr(0, inputfile.find_last_of('.'))
                : outputfile.substr(0, outputfile.find_last_of('.'))
            ) +  ".clone" + boost::lexical_cast<std::string>(rank) + ".h5";
            if (outputfile.empty())
                outputfile = inputfile.substr(0, inputfile.find_last_of('.')) +  ".out.h5";
        } catch (...) {
            std::cerr << "usage: [-T timelimit] [-c] inputfile [outputfile]" << std::endl;
            std::cerr << options << std::endl;
            std::abort();
        }
    }

    bool resume;
    std::size_t timelimit;
    std::string inputfile;
    std::string outputfile;
    std::string checkpointfile;
};

struct stop_callback {
    stop_callback(std::size_t timelimit)
        : limit(timelimit)
        , start(boost::chrono::high_resolution_clock::now())
    {}

    bool operator()() {
        return !signals.empty() 
            || (limit.count() > 0 && boost::chrono::high_resolution_clock::now() > start + limit);
    }

    boost::chrono::duration<std::size_t> limit;
    alps::ngs::signal signals;
    boost::chrono::high_resolution_clock::time_point start;
};

int main(int argc, char *argv[]) {

    try {
        boost::mpi::environment env(argc, argv);
        boost::mpi::communicator comm;

        args options(argc, argv, comm.rank());
        alps::parameters_type<ising_sim>::type parameters;
        if (comm.rank() == 0) {
            std::string suffix = options.inputfile.substr(options.inputfile.find_last_of('.'));
            if (suffix == ".xml")
                parameters = alps::make_parameters_from_xml(options.inputfile);
            else if (suffix == ".h5")
                alps::hdf5::archive(options.inputfile)["/parameters"] >> parameters;
            else
                throw std::runtime_error("Unsupported input formant: " + suffix + "!");
        }
        broadcast(comm, parameters);

        ising_sim sim(parameters, comm);

        if (options.resume && boost::filesystem::exists(options.checkpointfile))
            sim.load(options.checkpointfile);

        sim.run(stop_callback(options.timelimit));

        sim.save(options.checkpointfile);

        using alps::collect_results;
        alps::results_type<ising_sim>::type results = collect_results(sim);

        if (comm.rank() == 0) {
            std::cout << results << std::endl;
            alps::hdf5::archive ar(options.outputfile, "w");
            ar["/parameters"] << parameters;
            ar["/simulation/results"] << results;
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
