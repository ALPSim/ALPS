// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// pyalps.ngs.results: ALPSCore's result_set, as a mapping from names to
// pyalps.ngs.result (see accumulator_handles.hpp).
#define PY_ARRAY_UNIQUE_SYMBOL pyngsresults_PyArrayHandle
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <alps/hdf5.hpp>
#include "accumulator_handles.hpp"
#include "archive_savable.hpp"
#include <memory>
#include <sstream>
#include <string>
namespace nb = nanobind;
namespace acc = alps::accumulators;
namespace {
std::string print(acc::result_set const & self) {
    std::ostringstream out;
    out << self;
    return out.str();
}
// Loads the results at `path` into a new set first, so a failed load leaves
// the set and the archive's context as they were.
void load(acc::result_set & self, alps::hdf5::archive & ar, std::string const & path) {
    alps::hdf5::archive reader(ar);
    reader.set_context(ar.complete_path(path));
    acc::result_set loaded;
    loaded.load(reader);
    pyalps::replace(self, loaded);
}
} // namespace
NB_MODULE(pyngsresults_c, m) {
    nb::module_::import_("pyalps.cxx.pyngsresult_c");
    nb::class_<acc::result_set>(m, "results")
        .def(nb::init<>())
        .def("__len__",      [](acc::result_set const & self) { return self.size(); })
        .def("__contains__", [](acc::result_set const & self, std::string const & k) {
                                 return self.has(k);
                             })
        .def("__getitem__",  [](acc::result_set & self, std::string const & k) {
                                 return pyalps::result{pyalps::item(self, k)};
                             })
        .def("__setitem__",  [](acc::result_set & self, std::string const & k, pyalps::result const & v) {
                                 // The set and the assigned name share the result;
                                 // a result without measurements is refused.
                                 if (v.empty())
                                     throw nb::value_error("result has no measurements");
                                 pyalps::assign(self, k, v.value);
                             })
        .def("__delitem__",  &pyalps::erase<acc::result_set>)
        .def("__iter__",     [](acc::result_set & self) {
                                 nb::list keys;
                                 for (auto const & entry : self)
                                     keys.append(nb::cast(entry.first));
                                 return keys.attr("__iter__")();
                             })
        // keys/values/items are deliberately NOT defined here. Boost.Python's
        // map_indexing_suite did not define them either, so they resolved through
        // MutableMapping to set-like KeysView/ValuesView/ItemsView. Defining them
        // natively as nanobind iterators would narrow that surface (no len(), no
        // set operators, exhausted after one pass) and pyalps/ngs.py cannot
        // recover it -- its guard skips any name the C++ class already provides.
        .def("__str__",      &print)
        .def("save",         [](acc::result_set const & self, alps::hdf5::archive & ar) { self.save(ar); })
        .def("load",         &load);
    pyalps::mark_archive_savable(m.attr("results"));
}
