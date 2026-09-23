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
#include <nanobind/stl/string.h>
#include <alps/hdf5.hpp>
#include <alps/ngs/mcresults.hpp>
#include "../archive_savable.hpp"
#include "../mapping_lifetime.hpp"
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
            alps::hdf5::archive reader(ar);
            reader.set_context(ar.complete_path(path));
            alps::mcresults loaded;
            loaded.load(reader);
            while (!self.empty())
                pyalps::erase_map_item(self, self.begin()->first);
            self.swap(loaded);
        }
    }
}
NB_MODULE(pyngsresults_c, m) {
    nb::class_<alps::mcresults>(m, "results")
        .def(nb::init<>())
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
                                 if (self.has(k)) self[k] = v;
                                 else self.insert(k, v);
                             })
        .def("__delitem__", &pyalps::erase_map_item<alps::mcresults>)
        .def("__iter__",     [](alps::mcresults & self) {
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
        .def("__str__",      &alps::detail::mcresults_print)
        .def("save",         &alps::mcresults::save)
        .def("load",         &alps::detail::mcresults_load);
    pyalps::mark_archive_savable(m.attr("results"));
}
