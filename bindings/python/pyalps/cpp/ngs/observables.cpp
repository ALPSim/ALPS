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
// pyngsobservables_c — nanobind port.
//
// Unlike mcresults, alps::mcobservables doesn't override erase() so we
// can stand on a straight nb::class_<> with hand-written map methods
// that match the map_indexing_suite surface area. We keep it explicit
// (rather than nb::bind_map) because the class also carries non-map
// methods (reset/save/load/__lshift__/create*) that need to live on
// the same binding, and mixing bind_map with extra defs is noisy.
#define PY_ARRAY_UNIQUE_SYMBOL pyngsobservables_PyArrayHandle
#include <alps/hdf5/archive.hpp>
#include <alps/ngs/mcobservable.hpp>
#include <alps/ngs/mcobservables.hpp>
#include <alps/ngs/observablewrappers.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/make_iterator.h>
#include <nanobind/stl/string.h>
namespace nb = nanobind;
#include <cstdint>
#include <string>
namespace {
void mcobservables_load(alps::mcobservables & self, alps::hdf5::archive & ar, std::string const & path) {
    std::string current = ar.get_context();
    ar.set_context(path);
    self.load(ar);
    ar.set_context(current);
}
void createRealObservable(alps::mcobservables & self, std::string const & name, std::uint32_t binnum) {
    self << alps::ngs::RealObservable(name, binnum);
}
void createRealVectorObservable(alps::mcobservables & self, std::string const & name, std::uint32_t binnum) {
    self << alps::ngs::RealVectorObservable(name, binnum);
}
void addObservable(alps::mcobservables & self, nb::object const & obj) {
    // Mirror boost::python::call_method<void>(obj, "addToObservables", ref(self)):
    // bounce the call back into Python, passing `self` by reference.
    obj.attr("addToObservables")(nb::cast(&self, nb::rv_policy::reference));
}
}  // namespace
NB_MODULE(pyngsobservables_c, m) {
    nb::class_<alps::mcobservables>(m, "observables")
        .def(nb::init<>())
        .def("__len__",      [](alps::mcobservables const & self) { return self.size(); })
        .def("__contains__", [](alps::mcobservables const & self, std::string const & k) {
                                 return self.has(k);
                             })
        .def("__getitem__",  [](alps::mcobservables & self, std::string const & k) -> alps::mcobservable & {
                                 if (!self.has(k))
                                     throw nb::key_error(k.c_str());
                                 return self[k];
                             },
                             nb::rv_policy::reference_internal)
        .def("__setitem__",  [](alps::mcobservables & self, std::string const & k, alps::mcobservable const & v) {
                                 self.insert(k, v);
                             })
        .def("__iter__",     [](alps::mcobservables & self) {
                                 return nb::make_key_iterator(nb::type<alps::mcobservables>(), "key_iterator", self.begin(), self.end());
                             },
                             nb::keep_alive<0, 1>())
        .def("keys",         [](alps::mcobservables & self) {
                                 return nb::make_key_iterator(nb::type<alps::mcobservables>(), "key_iterator", self.begin(), self.end());
                             },
                             nb::keep_alive<0, 1>())
        .def("values",       [](alps::mcobservables & self) {
                                 return nb::make_value_iterator(nb::type<alps::mcobservables>(), "value_iterator", self.begin(), self.end());
                             },
                             nb::keep_alive<0, 1>())
        .def("items",        [](alps::mcobservables & self) {
                                 return nb::make_iterator(nb::type<alps::mcobservables>(), "item_iterator", self.begin(), self.end());
                             },
                             nb::keep_alive<0, 1>())
        .def("reset", &alps::mcobservables::reset, nb::arg("equilibrated") = false)
        .def("save",  &alps::mcobservables::save)
        .def("load",  &mcobservables_load)
        .def("__lshift__", &addObservable)
        .def("createRealObservable",       &createRealObservable,
             nb::arg("name"), nb::arg("binnum") = 0)
        .def("createRealVectorObservable", &createRealVectorObservable,
             nb::arg("name"), nb::arg("binnum") = 0);
}
