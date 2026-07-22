// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <nanobind/nanobind.h>
#include <nanobind/make_iterator.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <alps/hdf5/archive.hpp>
#include <alps/ngs/params.hpp>
#include <alps/ngs/detail/paramvalue.hpp>
#include <complex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <boost/variant/apply_visitor.hpp>
#include <boost/variant/static_visitor.hpp>
#include <vector>
namespace nb = nanobind;
namespace {
// Convert a Python dict into an alps::params. Same shape as the
// helper in mcbase.cpp but kept local to params.cpp so a change to
// the dispatch (e.g. adding complex support) can stay in one place
// alongside the other setitem logic.
alps::params py_dict_to_params(nb::dict const & d);
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
// via paramproxy's templated operator=.
void params_setitem(alps::params & self, nb::object const & key_obj, nb::object const & value) {
    std::string key = nb::cast<std::string>(nb::str(key_obj));
    if (nb::isinstance<nb::bool_>(value))
        self[key] = nb::cast<bool>(value);
    else if (nb::isinstance<nb::int_>(value))
        self[key] = nb::cast<int>(value);
    else if (nb::isinstance<nb::float_>(value))
        self[key] = nb::cast<double>(value);
    else if (nb::isinstance<nb::str>(value))
        self[key] = nb::cast<std::string>(value);
    else if (nb::isinstance<nb::list>(value) || nb::isinstance<nb::tuple>(value)) {
        // Heuristic: try doubles first, strings as fallback.
        try {
            self[key] = nb::cast<std::vector<double>>(value);
        } catch (nb::cast_error &) {
            self[key] = nb::cast<std::vector<std::string>>(value);
        }
    } else {
        throw nb::type_error("unsupported value type for params[]");
    }
}
nb::object params_getitem(alps::params & self, nb::object const & key_obj) {
    std::string key = nb::cast<std::string>(nb::str(key_obj));
    if (!self.defined(key))
        return nb::none();
    // params doesn't expose the underlying map directly, but
    // paramiterator yields (key, paramvalue) pairs; walk it to find the
    // entry and hand the variant to paramvalue_to_py.
    for (auto it = self.begin(); it != self.end(); ++it)
        if (it->first == key)
            return paramvalue_to_py(it->second);
    return nb::none();  // defensive — defined()==true should guarantee a hit
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
// Materialise an alps::params from a Python dict. Re-uses the same
// type dispatch as params_setitem so a round-tripped dict-built
// params contains exactly the same variant alternatives.
alps::params py_dict_to_params(nb::dict const & d) {
    alps::params p;
    for (auto item : d) {
        std::string k = nb::cast<std::string>(nb::str(item.first));
        nb::handle v = item.second;
        if (nb::isinstance<nb::bool_>(v))
            p[k] = nb::cast<bool>(v);
        else if (nb::isinstance<nb::int_>(v))
            p[k] = nb::cast<int>(v);
        else if (nb::isinstance<nb::float_>(v))
            p[k] = nb::cast<double>(v);
        else if (nb::isinstance<nb::str>(v))
            p[k] = nb::cast<std::string>(v);
        else if (nb::isinstance<nb::list>(v) || nb::isinstance<nb::tuple>(v)) {
            try { p[k] = nb::cast<std::vector<double>>(v); }
            catch (nb::cast_error &) {
                p[k] = nb::cast<std::vector<std::string>>(v);
            }
        } else {
            throw nb::type_error(
                ("unsupported value type for params key '" + k + "'").c_str());
        }
    }
    return p;
}
}  // namespace
NB_MODULE(pyngsparams_c, m) {
    nb::class_<alps::params>(m, "params")
        .def(nb::init<>())
        .def("__init__",
             [](alps::params * self, nb::dict const & d) {
                 new (self) alps::params(py_dict_to_params(d));
             },
             nb::arg("dict"))
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
