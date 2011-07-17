/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Matthias Troyer <troyer@comp-phys.org>             *
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

#define PY_ARRAY_UNIQUE_SYMBOL pyising_PyArrayHandle

#include "ising.hpp"

#include <alps/ngs/hdf5.hpp>
#include <alps/ngs/boost_python.hpp>

#include <alps/python/make_copy.hpp>

#include <boost/bind.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/return_internal_reference.hpp>

namespace alps {
    namespace detail {

        class ising_export : public ising_sim<alps::base> {

            public:

                ising_export(boost::python::object arg, std::size_t seed_offset = 42)
                    : ising_sim(alps::params(arg), seed_offset)
                {}

                bool run(boost::python::object stop_callback) {
                    ising_sim::run(boost::bind(ising_export::callback_wrapper, stop_callback));
                }

                typename base::parameters_type & get_params() {
                    return base::params;
                }

                mcobservables & get_measurements() {
                    return base::measurements;
                }

                double get_random() {
                    return alps::base::random();
                }

                void load(alps::hdf5::archive & ar, std::string const & path) {
                    std::string current = ar.get_context();
                    ar.set_context(path);
                    alps::base::load(ar);
                    ar.set_context(current);
                }

            private:

                static bool callback_wrapper(boost::python::object stop_callback) {
                   return boost::python::call<bool>(stop_callback.ptr());
                }

        };

    }
}

BOOST_PYTHON_MODULE(pyising_c) {

    boost::python::class_<alps::detail::ising_export, boost::noncopyable>(
          "sim",
          boost::python::init<boost::python::object, boost::python::optional<std::size_t> >()
    )
        .add_property("params", boost::python::make_function(&alps::detail::ising_export::get_params, boost::python::return_internal_reference<>()))
        .add_property("measurements", boost::python::make_function(&alps::detail::ising_export::get_measurements, boost::python::return_internal_reference<>()))
        .def("run", &alps::detail::ising_export::run)
        .def("random", &alps::detail::ising_export::get_random)
        .def("fraction_completed", boost::python::pure_virtual(&alps::detail::ising_export::fraction_completed))
        .def("save", &alps::detail::ising_export::save)
        .def("load", &alps::detail::ising_export::load)
    ;

}