// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
// <alps/ngs.hpp> umbrella retired in Phase 4 Slice 6 of the NGS
// retirement (ngs-retirement-scoping.md). This binding's surviving
// surface (`saveResults`) only needs `alps::mcresults`,
// `alps::params`, and `alps::hdf5::archive` — pull the narrow
// headers directly.
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <alps/ngs/mcresults.hpp>
#include <alps/ngs/params.hpp>
#include <alps/hdf5/archive.hpp>
#include <alps/mcbase.hpp>
namespace nb = nanobind;
namespace alps {
    namespace detail {
        void save_results_export(mcresults const & res, params const & par, alps::hdf5::archive & ar, std::string const & path) {
            ar["/parameters"] << par;
            if (res.size())
                ar[path] << res;
        }
    }
}
NB_MODULE(pyngsapi_c, m) {
    m.def("collectResults", [](alps::mcbase const & sim) {
        return alps::collect_results<alps::mcbase>(sim);
    });
    m.def("saveResults", &alps::detail::save_results_export);
}
