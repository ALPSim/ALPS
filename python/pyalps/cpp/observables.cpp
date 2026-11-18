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
// pyngsobservables_c: pyalps.ngs.observables, alps::alea::observable_set,
// as a mapping from names to pyalps.ngs.observable (see
// accumulator_handles.hpp).

#define PY_ARRAY_UNIQUE_SYMBOL pyngsobservables_PyArrayHandle

#include <alps/hdf5/archive.hpp>
#include "accumulator_handles.hpp"
#include "archive_savable.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
namespace nb = nanobind;

#include <cstdint>
#include <string>

namespace alea = alps::alea;

namespace {

void create(alea::observable_set & self, std::string const & name, std::uint32_t binnum) {
    self.insert(name, pyalps::real_observable(name, binnum));
}

void create_vector(alea::observable_set & self, std::string const & name, std::uint32_t binnum) {
    self.insert(name, pyalps::real_vector_observable(name, binnum));
}

void addObservable(alea::observable_set & self, nb::object const & obj) {
    // Mirror boost::python::call_method<void>(obj, "addToObservables", ref(self)):
    // bounce the call back into Python, passing `self` by reference.
    // nanobind's instance registry returns the already-registered
    // wrapper for &self (the one this call came through, verified
    // empirically), so the callback sees the identical Python object —
    // no duplicate wrapper, no separate lifetime to manage.
    obj.attr("addToObservables")(nb::cast(&self, nb::rv_policy::reference));
}

}  // namespace

NB_MODULE(pyngsobservables_c, m) {
    nb::module_::import_("pyalps.cxx.pyngsobservable_c");

    nb::class_<alea::observable_set>(m, "observables")
        .def(nb::init<>())
        .def("__len__",      [](alea::observable_set const & self) { return self.names().size(); })
        .def("__contains__", [](alea::observable_set const & self, std::string const & k) {
                                 return self.has(k);
                             })
        .def("__getitem__",  [](alea::observable_set & self, std::string const & k) {
                                 pyalps::need(self, k);
                                 return self[k];
                             })
        .def("__setitem__",  [](alea::observable_set & self, std::string const & k,
                                alea::accumulator const & v) { self.insert(k, v); })
        .def("__delitem__",  [](alea::observable_set & self, std::string const & k) {
                                 pyalps::need(self, k);
                                 self.erase(k);
                             })
        .def("__iter__",     &pyalps::names_iterator<alea::observable_set>)
        // keys/values/items are deliberately NOT defined here. Boost.Python's
        // map_indexing_suite did not define them either, so they resolved through
        // MutableMapping to set-like KeysView/ValuesView/ItemsView. Defining them
        // natively as nanobind iterators would narrow that surface (no len(), no
        // set operators, exhausted after one pass) and pyalps/ngs.py cannot
        // recover it -- its guard skips any name the C++ class already provides.
        // There is no equilibration flag; the argument is accepted for legacy
        // callers and ignored.
        .def("reset", [](alea::observable_set & self, bool /*equilibrated*/) { self.reset(); },
             nb::arg("equilibrated") = false)
        .def("save",  [](alea::observable_set const & self, alps::hdf5::archive & ar) { self.save(ar); })
        .def("load",  [](alea::observable_set & self, alps::hdf5::archive & ar, std::string const & path) {
                          self = pyalps::load_at<alea::observable_set>(ar, path);
                      })
        .def("__lshift__", &addObservable)
        .def("createRealObservable",       &create,
             nb::arg("name"), nb::arg("binnum") = 0)
        .def("createRealVectorObservable", &create_vector,
             nb::arg("name"), nb::arg("binnum") = 0);
    pyalps::mark_archive_savable(m.attr("observables"));
}
