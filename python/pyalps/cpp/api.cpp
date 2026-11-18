// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <alps/mc/api.hpp>
#include <alps/mc/mcbase.hpp>
#include <alps/params.hpp>
#include <alps/hdf5/archive.hpp>
namespace nb = nanobind;
NB_MODULE(pyngsapi_c, m) {
    nb::module_::import_("pyalps.cxx.pyngsresults_c");
    m.def("collectResults", [](alps::mcbase const & sim) { return alps::collect_results(sim); });
    m.def("saveResults", [](alps::alea::result_set const & res, alps::params const & par,
                            alps::hdf5::archive & ar, std::string const & path) {
        ar["/parameters"] << par;
        if (!res.names().empty())
            ar[path] << res;
    });
}
