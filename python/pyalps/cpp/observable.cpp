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
// pyngsobservable_c: pyalps.ngs.observable, an accumulator of ALPSCore's
// accumulators library (alps::accumulators::accumulator_wrapper).
#define PY_ARRAY_UNIQUE_SYMBOL pyngsobservable_PyArrayHandle
#include <alps/hdf5/archive.hpp>
#include "accumulator_handles.hpp"
#include "archive_savable.hpp"
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/stl/string.h>
namespace nb = nanobind;
#include <cstdint>
#include <memory>
#include <string>
#include <vector>
namespace acc = alps::accumulators;
namespace {
void observable_append(acc::accumulator_wrapper & self, nb::object const & data) {
    if (nb::isinstance<nb::int_>(data) || nb::isinstance<nb::float_>(data)) {
        self << nb::cast<double>(data);
        return;
    }
    try {
        auto values = nb::cast<nb::ndarray<double const, nb::c_contig>>(data);
        self << std::vector<double>(values.data(), values.data() + values.size());
    } catch (nb::cast_error const &) {
        throw nb::type_error("observable samples must be numeric scalars or contiguous float64 arrays");
    }
}
void observable_load(acc::accumulator_wrapper & self, alps::hdf5::archive & ar, std::string const & path) {
    std::string current = ar.get_context();
    ar.set_context(path);
    try {
        self.load(ar);
    } catch (...) {
        ar.set_context(current);
        throw;
    }
    ar.set_context(current);
}
} // namespace
NB_MODULE(pyngsobservable_c, m) {
    m.def("createRealObservable", &pyalps::full_binning<double>,
          nb::arg("name"), nb::arg("binnum") = 0);
    m.def("createRealVectorObservable", &pyalps::full_binning<std::vector<double>>,
          nb::arg("name"), nb::arg("binnum") = 0);
    nb::class_<acc::accumulator_wrapper>(m, "observable")
        .def("append",          &observable_append)
        // obs << value appends and returns obs so it chains. Bound in
        // C++ (rv_policy::none returns the existing wrapper) instead of
        // the former ngs.py monkeypatch onto the extension type, which
        // would break if nanobind ever marks its types immutable.
        .def("__lshift__",
             [](acc::accumulator_wrapper & self, nb::object const & data) -> acc::accumulator_wrapper & {
                 observable_append(self, data);
                 return self;
             },
             nb::rv_policy::none)
        .def_prop_ro("count", [](acc::accumulator_wrapper const & self) { return self.count(); })
        .def("merge",           [](acc::accumulator_wrapper & self, acc::accumulator_wrapper const & other) {
                                    self.merge(other);
                                })
        .def("save",            [](acc::accumulator_wrapper const & self, alps::hdf5::archive & ar) {
                                    self.save(ar);
                                })
        .def("load",            &observable_load)
        // Mirrors the legacy Boost.Python module, which (oddly, but
        // load-compatibly) bound addToObservable to the same helper
        // as load.
        .def("addToObservable", &observable_load);
    pyalps::mark_archive_savable(m.attr("observable"));
}
