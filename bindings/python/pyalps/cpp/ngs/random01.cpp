// Copyright (C) 2010 - 2013 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <alps/ngs/random01.hpp>
#include <alps/hdf5/archive.hpp>
namespace nb = nanobind;
NB_MODULE(pyngsrandom01_c, m) {
    nb::class_<alps::random01>(m, "random01")
        .def(nb::init<int>(), nb::arg("seed") = 42)
        .def("__deepcopy__",
             // copy.deepcopy() passes (self, memo); memo is unused.
             [](alps::random01 const & self, nb::handle /*memo*/) {
                 return alps::random01(self);
             })
        .def("__call__",
             static_cast<alps::random01::result_type (alps::random01::*)()>(
                 &alps::random01::operator()))
        .def("save", &alps::random01::save)
        .def("load", &alps::random01::load);
}
