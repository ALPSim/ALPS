/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Matthias Troyer <troyer@comp-phys.org>             *
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

#define PY_ARRAY_UNIQUE_SYMBOL pyngsparams_PyArrayHandle

#include <alps/hdf5/archive.hpp>
#include <alps/ngs/params.hpp>

#include <alps/python/make_copy.hpp>

#include <boost/python/tuple.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/to_python_converter.hpp>

#include <string>
#include <sstream>

namespace alps {
    namespace detail {

        std::size_t params_len(alps::params const & self) {
            return self.size();
        }

        boost::python::object params_getitem(alps::params & self, boost::python::object const & key) {
            if (self.defined(boost::python::call_method<std::string>(key.ptr(), "__str__")))
                return self[boost::python::call_method<std::string>(key.ptr(), "__str__")].cast<boost::python::object>();
            else
                return boost::python::object();
        }

        void params_setitem(alps::params & self, boost::python::object const & key, boost::python::object & value) {
            self[boost::python::call_method<std::string>(key.ptr(), "__str__")] = value;
        }

        void params_delitem(alps::params & self, boost::python::object const & key) {
            return self.erase(boost::python::call_method<std::string>(key.ptr(), "__str__"));
        }

        bool params_contains(alps::params & self, boost::python::object const & key) {
            return self.defined(boost::python::call_method<std::string>(key.ptr(), "__str__"));
        }

        boost::python::object value_or_default(alps::params & self, boost::python::object const & key, boost::python::object const & value) {
            return params_contains(self, key) ? params_getitem(self, key) : value;
        }

        void params_load(alps::params & self, alps::hdf5::archive & ar, std::string const & path = "/parameters") {
            std::string current = ar.get_context();
            ar.set_context(path);
            self.load(ar);
            ar.set_context(current);
        }
        BOOST_PYTHON_FUNCTION_OVERLOADS(params_load_overloads, params_load, 2, 3)

        struct param_iterator_to_python {
            static PyObject* convert(std::pair<std::string const, alps::detail::paramvalue> const & value) {
                return boost::python::incref(boost::python::str(value.first).ptr());
            }
        };

        boost::python::str params_print(alps::params & self) {
            std::stringstream ss;
            ss << self;
            return boost::python::str(ss.str());
        }

    }
}

BOOST_PYTHON_MODULE(pyngsparams_c) {

    boost::python::to_python_converter<
        std::pair<std::string const, alps::detail::paramvalue>,
        alps::detail::param_iterator_to_python
    >();

    boost::python::class_<alps::params>(
        "params",
        boost::python::init<boost::python::optional<boost::python::dict> >()
    )
        .def(boost::python::init<alps::hdf5::archive, boost::python::optional<std::string const &> >())
        .def(boost::python::init<boost::python::str const &>())

        .def("__len__", &alps::detail::params_len)
        .def("__deepcopy__", &alps::python::make_copy<alps::params>)
        .def("__getitem__", &alps::detail::params_getitem)
        .def("__setitem__", &alps::detail::params_setitem)
        .def("__delitem__", &alps::detail::params_delitem)
        .def("__contains__", &alps::detail::params_contains)
        .def("__iter__", boost::python::iterator<alps::params>())
        .def("__str__", &alps::detail::params_print)
        .def("valueOrDefault", &alps::detail::value_or_default)
        .def("save", &alps::params::save)
        .def("load", &alps::detail::params_load, alps::detail::params_load_overloads())
    ;
}
