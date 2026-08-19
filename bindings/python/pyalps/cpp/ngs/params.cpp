// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <nanobind/nanobind.h>
#include <nanobind/make_iterator.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <alps/hdf5/archive.hpp>
#include <alps/ngs/params.hpp>
#include <alps/ngs/detail/paramvalue.hpp>
#include <boost/filesystem/path.hpp>
#include <complex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/static_visitor.hpp>
#include <vector>
#include "../dict_to_params.hpp"
namespace nb = nanobind;
namespace {
// Walk the paramvalue variant and wrap each native alternative as a
// nb::object. Called from __getitem__.
struct paramvalue_to_py_visitor : boost::static_visitor<nb::object> {
    template <typename T>
    nb::object operator()(T const & value) const {
        return nb::cast(value);
    }
};
nb::object paramvalue_to_py(alps::detail::paramvalue const & pv) {
    return boost::apply_visitor(
        paramvalue_to_py_visitor(),
        static_cast<alps::detail::paramvalue_base const &>(pv));
}
// Deposit a native C++ value from a Python object into the paramvalue
// via paramproxy's templated operator= — shared ladder in
// ../dict_to_params.hpp so params, mcbase and the application modules
// all ingest values identically.
void params_setitem(alps::params & self, nb::object const & key_obj, nb::handle value) {
    pyalps::set_param_value(self, nb::cast<std::string>(nb::str(key_obj)), value);
}
nb::object params_getitem(alps::params & self, nb::object const & key_obj) {
    std::string key = nb::cast<std::string>(nb::str(key_obj));
    alps::detail::paramvalue const * value = self.find(key);
    return value ? paramvalue_to_py(*value) : nb::none();
}
void params_delitem(alps::params & self, nb::object const & key_obj) {
    self.erase(nb::cast<std::string>(nb::str(key_obj)));
}
bool params_contains(alps::params & self, nb::object const & key_obj) {
    return self.defined(nb::cast<std::string>(nb::str(key_obj)));
}
nb::object value_or_default(alps::params & self, nb::object const & key, nb::object const & dflt) {
    return params_contains(self, key) ? params_getitem(self, key) : dflt;
}
void params_load(alps::params & self, alps::hdf5::archive & ar, std::string const & path) {
    std::string current = ar.get_context();
    ar.set_context(path);
    self.load(ar);
    ar.set_context(current);
}
std::string params_print(alps::params & self) {
    std::stringstream ss;
    ss << self;
    return ss.str();
}
// deepcopy support — nanobind passes (self, memo); memo unused.
alps::params params_deepcopy(alps::params const & self, nb::handle /*memo*/) {
    return alps::params(self);
}
}  // namespace
NB_MODULE(pyngsparams_c, m) {
    nb::class_<alps::params>(m, "params")
        .def(nb::init<>())
        .def("__init__",
             [](alps::params * self, nb::dict const & d) {
                 new (self) alps::params(pyalps::params_from_dict(d));
             },
             nb::arg("dict"))
        // Read a classic ALPS text parameter file, matching the str
        // constructor of the Boost.Python module.
        .def("__init__",
             [](alps::params * self, std::string const & filename) {
                 new (self) alps::params(boost::filesystem::path(filename));
             },
             nb::arg("filename"))
        .def(nb::init<alps::hdf5::archive, std::string const &>(),
             nb::arg("archive"),
             nb::arg("path") = std::string("/parameters"))
        .def("__len__",      [](alps::params const & self) { return self.size(); })
        .def("__deepcopy__", &params_deepcopy)
        .def("__getitem__",  &params_getitem)
        .def("__setitem__",  &params_setitem)
        .def("__delitem__",  &params_delitem)
        .def("__contains__", &params_contains)
        .def("__iter__",     [](alps::params & self) {
                                 // paramiterator yields pair<string const, paramvalue>;
                                 // make_key_iterator projects out pair.first.
                                 return nb::make_key_iterator(
                                     nb::type<alps::params>(),
                                     "key_iterator",
                                     self.begin(), self.end());
                             },
                             nb::keep_alive<0, 1>())
        .def("__str__",      &params_print)
        .def("valueOrDefault", &value_or_default)
        .def("save",         &alps::params::save)
        .def("load",         &params_load,
             nb::arg("archive"),
             nb::arg("path") = std::string("/parameters"));
}
