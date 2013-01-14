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

#include "ising.hpp"

#include <alps/python/make_copy.hpp>
#include <alps/random/mersenne_twister.hpp> // TODO: why do we need this?

#include <boost/bind.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>

class ising_export : public ising_sim {

    public:

        ising_export(parameters_type const & params)
            : ising_sim(params)
        {}

        results_type collect_results(result_names_type const & names = result_names_type()) {
            return names.size() ? ising_sim::collect_results(names) : ising_sim::collect_results();
        }

        bool run(boost::python::object stop_callback) {
            return ising_sim::run(boost::bind(&ising_export::run_helper, this, stop_callback));
        }

        double get_random() {
            return ising_sim::random();
        }

        parameters_type & get_parameters() {
            return ising_sim::parameters;
        }

    private:

        bool run_helper(boost::python::object stop_callback) {
            return boost::python::call<bool>(stop_callback.ptr());
        }
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(collect_results_overloads, collect_results, 0, 1)

BOOST_PYTHON_MODULE(ising_c) {
    boost::python::class_<ising_export>(
          "sim",
          boost::python::init<ising_export::parameters_type const &>()
    )
        .add_property("random", &ising_export::get_random)
        .add_property("parameters", boost::python::make_function(&ising_export::get_parameters, boost::python::return_internal_reference<>()))
        .def("fraction_completed", &ising_export::fraction_completed)
        .def("run", &ising_export::run)
        .def("resultNames", &ising_export::result_names)
        .def("unsavedResultNames", &ising_export::unsaved_result_names)
        .def("collectResults", &ising_export::collect_results, collect_results_overloads(boost::python::args("names")))
//        .def("save", &ising_export::save) // TODO: implement
//        .def("load", &ising_export::load) // TODO: implement
    ;
}
