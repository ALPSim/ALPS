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

#define PY_ARRAY_UNIQUE_SYMBOL isingsim_PyArrayHandle

#include <alps/ngs/api.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/ngs/config.hpp>
#include <alps/ngs/signal.hpp>
#include <alps/ngs/params.hpp>
#include <alps/ngs/mcresults.hpp>
#include <alps/ngs/boost_python.hpp>
#include <alps/ngs/mcobservables.hpp>
#include <alps/ngs/observablewrappers.hpp>
#include <alps/ngs/make_parameters_from_xml.hpp>

#include <alps/python/make_copy.hpp>
#include <alps/random/mersenne_twister.hpp> // TODO: why do we need this?

#include <boost/bind.hpp>
#include <boost/chrono.hpp>
#include <boost/function.hpp>
#include <boost/python/dict.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/python/overloads.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/python/return_internal_reference.hpp>

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

        ising_sim(parameters_type const & parameters)
            : params(parameters)
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

        bool run(boost::function<bool ()> const & stop_callback) {
          bool stopped = false;
          do {
                update();
                measure();
            } while(!(stopped = stop_callback()) && fraction_completed() < 1.);
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

        parameters_type const & parameters() const { return params; }

        void save(alps::hdf5::archive & ar) const {
            ar["/parameters"] << params;
            std::string context = ar.get_context();
            ar.set_context("/simulation/realizations/0/clones/0");
            ar["measurements"] << measurements;

            ar.set_context("checkpoint");
            ar["length"] << length;
            ar["sweeps"] << sweeps;
            ar["thermalization_sweeps"] << thermalization_sweeps;
            ar["beta"] << beta;
            ar["spins"] << spins;
            ar["measurements"] << measurements;

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
            ar.set_context("/simulation/realizations/0/clones/0");
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

        friend double py_random(ising_sim &);
        friend parameters_type py_parameters(ising_sim &);
        friend observables_type py_measurements(ising_sim &);

        parameters_type params;
        boost::variate_generator<boost::mt19937, boost::uniform_real<> > mutable random;
        observables_type measurements;

        int length;
        int sweeps;
        int thermalization_sweeps;
        int total_sweeps;
        double beta;
        std::vector<int> spins;
};

bool py_run_helper(boost::python::object stop_callback) {
    return boost::python::call<bool>(stop_callback.ptr());
}
bool py_run(ising_sim & self, boost::python::object stop_callback) {
    return self.run(boost::bind(py_run_helper, stop_callback));
}

ising_sim::results_type py_collect_results(ising_sim & self, ising_sim::result_names_type const & names = ising_sim::result_names_type()) {
    return names.size() ? self.collect_results(names) : self.collect_results();
}
BOOST_PYTHON_FUNCTION_OVERLOADS(py_collect_results_overloads, py_collect_results, 1, 2)

double py_random(ising_sim & self) {
    return self.random();
}

void py_save(ising_sim const & self, alps::hdf5::archive & ar) {
    self.save(ar);
}

void py_load(ising_sim & self, alps::hdf5::archive & ar) {
    self.load(ar);
}

// TODO: make this nicer ...
BOOST_PYTHON_MODULE(exampleising_c) {
    boost::python::class_<ising_sim>(
          "sim",
          boost::python::init<ising_sim::parameters_type const &>()
    )
        //jan: why not restrict exports to necessary things? parameters, results, run, load, save?
        .add_property("random", &py_random)
        .add_property("parameters", boost::python::make_function(&ising_sim::parameters, boost::python::return_internal_reference<>()))
        .def("fraction_completed", &ising_sim::fraction_completed)
        .def("run", &py_run)
        .def("resultNames", &ising_sim::result_names)
        .def("unsavedResultNames", &ising_sim::unsaved_result_names)
        .def("collectResults", &py_collect_results, py_collect_results_overloads(boost::python::args("names")))
        .def("save", &py_save)
        .def("load", &py_load)
    ;
}
