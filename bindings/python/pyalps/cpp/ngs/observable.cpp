/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Matthias Troyer <troyer@comp-phys.org>             *
 *               2026       by the ALPS collaboration                              *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the "Software"),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
// pyngsobservable_c — nanobind port.
#define PY_ARRAY_UNIQUE_SYMBOL pyngsobservable_PyArrayHandle
#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/complex.hpp>
#include <alps/ngs/mcobservable.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
namespace nb = nanobind;
#include <alps/alea/detailedbinning.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <valarray>
namespace alps {
    namespace detail {
        void observable_append(alps::mcobservable & self, nb::object const & data) {
            if (nb::isinstance<nb::int_>(data) || nb::isinstance<nb::float_>(data)) {
                self << nb::cast<double>(data);
                return;
            }
            try {
                auto values = nb::cast<nb::ndarray<double, nb::c_contig>>(data);
                self << std::valarray<double>(values.data(), values.size());
            } catch (nb::cast_error const &) {
                throw nb::type_error("observable samples must be numeric scalars or contiguous float64 arrays");
            }
        }
        void observable_load(alps::mcobservable & self, alps::hdf5::archive & ar, std::string const & path) {
            std::string current = ar.get_context();
            ar.set_context(path);
            self.load(ar);
            ar.set_context(current);
        }
        alps::mcobservable create_RealObservable_export(std::string name) {
            return alps::mcobservable(std::make_shared<alps::RealObservable>(name).get());
        }
        alps::mcobservable create_RealVectorObservable_export(std::string name) {
            return alps::mcobservable(std::make_shared<alps::RealVectorObservable>(name).get());
        }
    }
}
NB_MODULE(pyngsobservable_c, m) {
    m.def("createRealObservable",       &alps::detail::create_RealObservable_export);
    m.def("createRealVectorObservable", &alps::detail::create_RealVectorObservable_export);
    nb::class_<alps::mcobservable>(m, "observable")
        .def("append",          &alps::detail::observable_append)
        .def("merge",           &alps::mcobservable::merge)
        .def("save",            &alps::mcobservable::save)
        .def("load",            &alps::detail::observable_load)
        .def("addToObservable", &alps::detail::observable_load);
}
