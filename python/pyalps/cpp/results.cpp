// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// pyalps.ngs.results: alps::alea::result_set, as a mapping from names to
// pyalps.ngs.result (see accumulator_handles.hpp).

#define PY_ARRAY_UNIQUE_SYMBOL pyngsresults_PyArrayHandle

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <alps/hdf5.hpp>
#include "accumulator_handles.hpp"
#include "archive_savable.hpp"

#include <sstream>
#include <string>

namespace nb = nanobind;
namespace alea = alps::alea;

namespace {

std::string print(alea::result_set const & self) {
    std::ostringstream out;
    out << self;
    return out.str();
}

} // namespace

NB_MODULE(pyngsresults_c, m) {
    nb::module_::import_("pyalps.cxx.pyngsresult_c");

    nb::class_<alea::result_set>(m, "results")
        .def(nb::init<>())
        .def("__len__",      [](alea::result_set const & self) { return self.names().size(); })
        .def("__contains__", [](alea::result_set const & self, std::string const & k) {
                                 return self.has(k);
                             })
        .def("__getitem__",  [](alea::result_set const & self, std::string const & k) {
                                 pyalps::need(self, k);
                                 return pyalps::wrap(self[k]);
                             })
        .def("__setitem__",  [](alea::result_set & self, std::string const & k, pyalps::result const & v) {
                                 // A result without measurements is refused.
                                 if (v.empty())
                                     throw nb::value_error("result has no measurements");
                                 self.insert(k, *v.value);
                             })
        .def("__delitem__",  [](alea::result_set & self, std::string const & k) {
                                 pyalps::need(self, k);
                                 self.erase(k);
                             })
        .def("__iter__",     &pyalps::names_iterator<alea::result_set>)
        // keys/values/items are deliberately NOT defined here. Boost.Python's
        // map_indexing_suite did not define them either, so they resolved through
        // MutableMapping to set-like KeysView/ValuesView/ItemsView. Defining them
        // natively as nanobind iterators would narrow that surface (no len(), no
        // set operators, exhausted after one pass) and pyalps/ngs.py cannot
        // recover it -- its guard skips any name the C++ class already provides.
        .def("__str__",      &print)
        .def("save",         [](alea::result_set const & self, alps::hdf5::archive & ar) { self.save(ar); })
        .def("load",         [](alea::result_set & self, alps::hdf5::archive & ar, std::string const & path) {
                                 self = pyalps::load_at<alea::result_set>(ar, path);
                             });
    pyalps::mark_archive_savable(m.attr("results"));
}
