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

#define PY_ARRAY_UNIQUE_SYMBOL pyngsresults_PyArrayHandle

#include <alps/hdf5.hpp>
#include <alps/ngs/mcresults.hpp>

#include <alps/ngs/boost_python.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>

namespace alps {
    namespace detail {

        std::string mcresults_print(alps::mcresults & self) {
            std::stringstream sstr;
            sstr << self;
            return sstr.str();
        }

        void mcresults_load(alps::mcresults & self, alps::hdf5::archive & ar, std::string const & path) {
            std::string current = ar.get_context();
            ar.set_context(path);
            self.load(ar);
            ar.set_context(current);
        }
    }
}

BOOST_PYTHON_MODULE(pyngsresults_c) {
    boost::python::class_<alps::mcresults>(
        "results",
        boost::python::no_init
    )
        .def(boost::python::map_indexing_suite<alps::mcresults>())
        .def("__str__", &alps::detail::mcresults_print)
        .def("save", &alps::mcresults::save)
        .def("load", &alps::detail::mcresults_load)
    ;

}
