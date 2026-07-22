/*****************************************************************************
*
* ALPS Project Applications: Directed Worm Algorithm
*
* Copyright (C) 2013 by Matthias Troyer  <troyer@phys.ethz.ch> ,
*                       Lode Pollet      <pollet@phys.ethz.ch> ,
*                       Ping Nang Ma     <pingnang@phys.ethz.ch>
*               2026 by the ALPS collaboration
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/
// dwa_c — nanobind port.
//
// The old boost::python binding registered each std::vector<T> we
// handed out as a distinct nb::class_<std::vector<T>> carrying the
// vector_indexing_suite. nanobind's built-in STL caster auto-
// converts std::vector<T> ↔ Python list for us, so those
// registrations are retired; any code that wrote
// `dwa_c.std_vector_double(...)` now just passes / receives a Python
// list directly.
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
namespace nb = nanobind;
#include "worldlines.hpp"
#include "bandstructure.hpp"
#include <string>
#include <vector>
NB_MODULE(dwa_c, m) {
    m.doc() = "ALPS DWA (directed worm algorithm) Python bindings.";
    nb::class_<kink>(m, "kink")
        .def(nb::init<unsigned int>(), nb::arg("siteindicator"))
        .def(nb::init<unsigned int, double, unsigned short>(),
             nb::arg("siteindicator"), nb::arg("time"), nb::arg("state"))
        .def("__repr__",     &kink::representation)
        .def("siteindicator",&kink::siteindicator)
        .def("time",         &kink::time)
        .def("state",        &kink::state);
    nb::class_<worldlines::location_type>(m, "location_type");
    nb::class_<worldlines>(m, "worldlines")
        .def(nb::init<>())
        .def(nb::init<unsigned int>(), nb::arg("num_sites"))
        .def("__repr__", &worldlines::representation)
        .def("load", static_cast<void (worldlines::*)(std::string const &)>(&worldlines::load))
        .def("save", static_cast<void (worldlines::*)(std::string const &) const>(&worldlines::save))
        .def("open_worldlines", &worldlines::open_worldlines)
        .def("worldlines_siteindicator", &worldlines::worldlines_siteindicator)
        .def("worldlines_time",          &worldlines::worldlines_time)
        .def("worldlines_state",         &worldlines::worldlines_state)
        .def("num_sites", &worldlines::num_sites)
        .def("num_kinks", &worldlines::num_kinks)
        .def("states",   &worldlines::states)
        .def("location", &worldlines::location)
        .def("state_before", &worldlines::state_before)
        .def("state",        &worldlines::state)
        .def("is_valid", static_cast<bool (worldlines::*)(unsigned short)>(&worldlines::is_valid));
    nb::class_<wormpair>(m, "wormpair")
        .def(nb::init<>())
        .def(nb::init<worldlines::location_type, kink, bool, bool>())
        .def("__repr__", &wormpair::representation)
        .def("wormhead", &wormpair::wormhead)
        .def("wormtail", &wormpair::wormtail)
        .def("wormhead_site",    &wormpair::site)
        .def("wormhead_time",    &wormpair::time)
        .def("wormhead_forward", &wormpair::forward)
        .def("wormtail_site", &wormpair::wormtail_site)
        .def("wormtail_time", &wormpair::wormtail_time)
        .def("next_partnersite", &wormpair::next_partnersite)
        .def("next_time",        &wormpair::next_time)
        .def("wormhead_turns_around",                         &wormpair::wormhead_turns_around)
        .def("wormhead_moves_to_new_time",                    &wormpair::wormhead_moves_to_new_time)
        .def("wormhead_inserts_vertex_and_jumps_to_new_site", &wormpair::wormhead_inserts_vertex_and_jumps_to_new_site)
        .def("wormhead_deletes_vertex_and_jumps_to_new_site", &wormpair::wormhead_deletes_vertex_and_jumps_to_new_site)
        .def("wormhead_relinks_vertex_and_jumps_to_new_site", &wormpair::wormhead_relinks_vertex_and_jumps_to_new_site)
        .def("wormhead_crosses_vertex",                       &wormpair::wormhead_crosses_vertex)
        .def("wormhead_annihilates_wormtail",                 &wormpair::wormhead_annihilates_wormtail);
    nb::class_<bandstructure>(m, "bandstructure")
        .def(nb::init<double, double, double, double, unsigned int>(),
             nb::arg("V0"), nb::arg("lambda"), nb::arg("a"), nb::arg("m"), nb::arg("L"))
        .def(nb::init<std::vector<double> const &, std::vector<double> const &,
                      double, double, unsigned int>(),
             nb::arg("V0"), nb::arg("lambda"), nb::arg("a"), nb::arg("m"), nb::arg("L"))
        .def("__repr__", static_cast<std::string (bandstructure::*)()>(&bandstructure::representation))
        .def("t",       static_cast<std::vector<double> (bandstructure::*)()>(&bandstructure::get_t))
        .def("U",       static_cast<double              (bandstructure::*)()>(&bandstructure::get_U))
        .def("Ut",      static_cast<std::vector<double> (bandstructure::*)()>(&bandstructure::get_Ut))
        .def("norm",    static_cast<std::vector<double> (bandstructure::*)()>(&bandstructure::get_norm))
        .def("q",       static_cast<std::vector<double> (bandstructure::*)(unsigned int)>(&bandstructure::get_q))
        .def("wk2",     static_cast<std::vector<double> (bandstructure::*)(unsigned int)>(&bandstructure::get_wk2))
        .def("wk2_c",   static_cast<double              (bandstructure::*)()>(&bandstructure::get_wk2_c))
        .def("wk2_d",   static_cast<double              (bandstructure::*)()>(&bandstructure::get_wk2_d));
}
