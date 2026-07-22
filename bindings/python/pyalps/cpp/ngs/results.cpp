// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
// A generic map binder cannot be used for alps::mcresults because
// mcresults::erase(std::string const &) shadows the std::map::erase(iterator)
// that bind_map relies on for __delitem__. Synthesise the dict-like surface
// by hand instead. (Same applies to nanobind's bind_map.)
#define PY_ARRAY_UNIQUE_SYMBOL pyngsresults_PyArrayHandle
#include <nanobind/nanobind.h>
#include <nanobind/make_iterator.h>
#include <nanobind/stl/string.h>
#include <alps/hdf5.hpp>
#include <alps/ngs/mcresults.hpp>
#include <sstream>
#include <stdexcept>
#include <string>
namespace nb = nanobind;
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
NB_MODULE(pyngsresults_c, m) {
    nb::class_<alps::mcresults>(m, "results")
        .def("__len__",      [](alps::mcresults const & self) { return self.size(); })
        .def("__contains__", [](alps::mcresults const & self, std::string const & k) {
                                 return self.has(k);
                             })
        .def("__getitem__",  [](alps::mcresults & self, std::string const & k) -> alps::mcresult const & {
                                 if (!self.has(k))
                                     throw nb::key_error(k.c_str());
                                 return self[k];
                             },
                             nb::rv_policy::reference_internal)
        .def("__setitem__",  [](alps::mcresults & self, std::string const & k, alps::mcresult const & v) {
                                 self.insert(k, v);
                             })
        .def("__delitem__",  [](alps::mcresults & self, std::string const & k) {
                                 if (!self.has(k))
                                     throw nb::key_error(k.c_str());
                                 self.erase(k);
                             })
        .def("__iter__",     [](alps::mcresults & self) {
                                 return nb::make_key_iterator(
                                     nb::type<alps::mcresults>(),
                                     "key_iterator",
                                     self.begin(), self.end());
                             },
                             nb::keep_alive<0, 1>())
        .def("keys",         [](alps::mcresults & self) {
                                 return nb::make_key_iterator(
                                     nb::type<alps::mcresults>(),
                                     "key_iterator",
                                     self.begin(), self.end());
                             },
                             nb::keep_alive<0, 1>())
        .def("values",       [](alps::mcresults & self) {
                                 return nb::make_value_iterator(
                                     nb::type<alps::mcresults>(),
                                     "value_iterator",
                                     self.begin(), self.end());
                             },
                             nb::keep_alive<0, 1>())
        .def("items",        [](alps::mcresults & self) {
                                 return nb::make_iterator(
                                     nb::type<alps::mcresults>(),
                                     "item_iterator",
                                     self.begin(), self.end());
                             },
                             nb::keep_alive<0, 1>())
        .def("__str__",      &alps::detail::mcresults_print)
        .def("save",         &alps::mcresults::save)
        .def("load",         &alps::detail::mcresults_load);
}
