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
// pyngsobservable_c: pyalps.ngs.observable, an observable of
// <alps/alea/observables.hpp> (alps::alea::accumulator).

#define PY_ARRAY_UNIQUE_SYMBOL pyngsobservable_PyArrayHandle

#include <alps/hdf5/archive.hpp>
#include "accumulator_handles.hpp"
#include "archive_savable.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
namespace nb = nanobind;

#include <cstdint>
#include <exception>
#include <string>
#include <vector>

namespace alea = alps::alea;

namespace {

// Appends a sample. One that does not fit the observable, such as a vector of
// another size, is a ValueError.
void append_sample(alea::accumulator & self, double const * values, std::size_t size) {
    try {
        self.add(values, size);
    } catch (std::exception const & error) {
        throw nb::value_error(error.what());
    }
}

void observable_append(alea::accumulator & self, nb::object const & data) {
    if (nb::isinstance<nb::int_>(data) || nb::isinstance<nb::float_>(data)) {
        double const value = nb::cast<double>(data);
        append_sample(self, &value, 1);
        return;
    }
    try {
        auto array = nb::cast<nb::ndarray<double const, nb::c_contig>>(data);
        append_sample(self, array.data(), array.size());
    } catch (nb::cast_error const &) {
        throw nb::type_error("observable samples must be numeric scalars or contiguous float64 arrays");
    }
}

// Replaces the measurements of `self` by those stored at `path`, for every
// reference to the observable.
void observable_load(alea::accumulator & self, alps::hdf5::archive & ar, std::string const & path) {
    alea::observable_set const loaded = pyalps::load_at<alea::observable_set>(ar, path);
    pyalps::need(loaded, self.name());
    alea::observable_set mine = pyalps::alone(self);
    self.reset();
    mine.merge(loaded.results().select({self.name()}));
}

} // namespace

NB_MODULE(pyngsobservable_c, m) {
    m.def("createRealObservable", &pyalps::real_observable,
          nb::arg("name"), nb::arg("binnum") = 0);
    m.def("createRealVectorObservable", &pyalps::real_vector_observable,
          nb::arg("name"), nb::arg("binnum") = 0);

    nb::class_<alea::accumulator>(m, "observable")
        .def("append",          &observable_append)
        // obs << value appends and returns obs so it chains. Bound in
        // C++ (rv_policy::none returns the existing wrapper) instead of
        // the former ngs.py monkeypatch onto the extension type, which
        // would break if nanobind ever marks its types immutable.
        .def("__lshift__",
             [](alea::accumulator & self, nb::object const & data) -> alea::accumulator & {
                 observable_append(self, data);
                 return self;
             },
             nb::rv_policy::none)
        .def_prop_ro("count", [](alea::accumulator const & self) { return self.count(); })
        // Appends the measurements of an independent run.
        .def("merge",           [](alea::accumulator & self, alea::accumulator const & other) {
                                    alea::observable_set mine = pyalps::alone(self), theirs;
                                    theirs.insert(self.name(), other);
                                    mine.merge(theirs.results());
                                })
        // The state to continue from, as a set of this one observable.
        .def("save",            [](alea::accumulator const & self, alps::hdf5::archive & ar) {
                                    pyalps::alone(self).save(ar);
                                })
        .def("load",            &observable_load)
        // Mirrors the legacy Boost.Python module, which (oddly, but
        // load-compatibly) bound addToObservable to the same helper
        // as load.
        .def("addToObservable", &observable_load);
    pyalps::mark_archive_savable(m.attr("observable"));
}
