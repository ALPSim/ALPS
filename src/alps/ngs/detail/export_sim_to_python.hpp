/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_DETAIL_EXPORT_SIM_TO_PYTHON_HPP
#define ALPS_NGS_DETAIL_EXPORT_SIM_TO_PYTHON_HPP

#include <alps/ngs/boost_python.hpp>

#include <alps/python/make_copy.hpp>

#include <boost/bind.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/wrapper.hpp>
#include <boost/python/return_internal_reference.hpp>

namespace alps {

    template<typename T> class export2python_wrapper : public T {

        public:

            export2python_wrapper(typename T::parameters_type const & parm, std::size_t seed_offset = 0)
                : T(parm, seed_offset)
            {}

            typename T::results_type collect_results(typename T::result_names_type const & names = typename T::result_names_type()) {
                return names.size() ? T::collect_results(names) : T::collect_results();
            }

            bool run(boost::python::object stop_callback) {
                return T::run(boost::bind(&export2python_wrapper<T>::run_helper, this, stop_callback));
            }

            alps::random01 & get_random() {
                return T::random;
            }

            typename T::parameters_type & get_parameters() {
                return T::parameters;
            }

        private:

            bool run_helper(boost::python::object stop_callback) {
              return boost::python::call<bool>(stop_callback.ptr());
            }
    };

}
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(collect_results_overloads, collect_results, 0, 1)

#define ALPS_EXPORT_SIM_TO_PYTHON(NAME, CLASS)                                                                                                              \
    boost::python::class_< alps::export2python_wrapper< CLASS >, boost::noncopyable, boost::python::bases<alps::mcbase> >(                                  \
          #NAME ,                                                                                                                                           \
          boost::python::init< CLASS ::parameters_type const &, boost::python::optional<std::size_t> >()                                                    \
    )                                                                                                                                                       \
        .add_property("random", boost::python::make_function(                                                                                               \
            &alps::export2python_wrapper< CLASS >::get_random, boost::python::return_internal_reference<>())                                                \
         )                                                                                                                                                  \
        .add_property("parameters", boost::python::make_function(                                                                                           \
            &alps::export2python_wrapper< CLASS >::get_parameters, boost::python::return_internal_reference<>())                                            \
         )                                                                                                                                                  \
        .def("run", static_cast<bool(alps::export2python_wrapper< CLASS >::*)(boost::python::object)>(&alps::export2python_wrapper< CLASS >::run))          \
        .def("resultNames", &alps::export2python_wrapper< CLASS >::result_names)                                                                            \
        .def("unsavedResultNames", &alps::export2python_wrapper< CLASS >::unsaved_result_names)                                                             \
        .def("collectResults", &alps::export2python_wrapper< CLASS >::collect_results, collect_results_overloads(boost::python::args("names")))             \
        .def("save", static_cast<void(alps::export2python_wrapper< CLASS >::*)(alps::hdf5::archive &) const>(&alps::export2python_wrapper< CLASS >::save))  \
        .def("load", static_cast<void(alps::export2python_wrapper< CLASS >::*)(alps::hdf5::archive &)>(&alps::export2python_wrapper< CLASS >::load))

#endif
