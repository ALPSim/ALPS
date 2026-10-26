// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <alps/hdf5/archive.hpp>
#include "archive_savable.hpp"
#include <alps/ngs/make_parameters_from_file.hpp>
#include <alps/params.hpp>
#include <boost/variant/static_visitor.hpp>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include "dict_to_params.hpp"
namespace nb = nanobind;
// pyalps.ngs.params wraps ALPSCore's alps::params. Values from Python are
// supplied as text (dict_to_params.hpp). A parameter the program has defined
// reads back with its defined type; one that is only supplied reads back
// from its text as params::save() writes it: an integer, a float, an array
// of numbers, or else the text (alps::params_ns::supplied_value).
namespace {
// The NumPy dtype of a vector element; long and unsigned long are not the
// fixed-width aliases on every platform, so the name comes from the size.
template <typename T>
std::string element_dtype() {
    if constexpr (std::is_same<T, bool>::value)
        return "bool";
    else if constexpr (std::is_floating_point<T>::value)
        return alps::python::numpy_dtype<T>::name;
    else
        return (std::is_signed<T>::value ? "int" : "uint") + std::to_string(8 * sizeof(T));
}
struct dict_value_to_py : boost::static_visitor<nb::object> {
    nb::object operator()(alps::params_ns::detail::None const &) const {
        return nb::none();
    }
    template <typename T>
    nb::object operator()(T const & value) const {
        return nb::cast(value);
    }
    template <typename T>
    nb::object operator()(std::vector<T> const & value) const {
        // String parameters need variable-length elements: NumPy's inferred
        // fixed-width Unicode dtype silently truncates longer replacements.
        if constexpr (std::is_same<T, std::string>::value)
            return nb::cast(value);
        // An empty numeric sequence has no elements from which NumPy can
        // infer its type. Preserve the native family, including Boolean masks.
        else
            return alps::python::numpy_module().attr("array")(
                nb::cast(value), nb::arg("dtype") = element_dtype<T>());
    }
};
bool params_contains(alps::params const & self, std::string const & key) {
    return self.exists(key) || self.supplied(key);
}
nb::object params_getitem(alps::params const & self, nb::object const & key_obj) {
    std::string const key = nb::cast<std::string>(nb::str(key_obj));
    if (self.exists(key))
        return self[key].apply_visitor(dict_value_to_py());
    auto const supplied = self.supplied_text().find(key);
    if (supplied != self.supplied_text().end())
        return alps::params_ns::supplied_value(key, supplied->second).apply_visitor(dict_value_to_py());
    throw nb::key_error(key.c_str());
}
// A copy of `self` without `key`, which params cannot remove. Its other
// values are kept; parameters lose their definitions and supplied text
// stays supplied, so programs define them again.
alps::params without(alps::params const & self, std::string const & key) {
    alps::params result;
    for (auto const & supplied : self.supplied_text())
        if (supplied.first != key)
            result.supply(supplied.first, supplied.second);
    for (auto const & entry : self)
        if (entry.first != key && !self.supplied(entry.first))
            result[entry.first] = entry.second;
    return result;
}
// The text parsed with the type of the value a parameter holds, as
// define<T>() would parse it; false for a parameter without a value.
struct parse_as : boost::static_visitor<bool> {
    alps::params & parsed;
    std::string const & key;
    bool operator()(alps::params_ns::detail::None const &) const {
        return false;
    }
    template <typename T>
    bool operator()(T const &) const {
        parsed.define<T>(key, "");
        if (parsed.has_missing())
            throw nb::type_error(("parameter '" + key + "' holds a "
                + alps::params_ns::detail::type_info<T>::pretty_name()
                + "; cannot assign '" + parsed.supplied_text().at(key) + "'").c_str());
        return true;
    }
};
// A parameter that holds a value takes the new one with its type; any other
// is supplied as text.
void params_setitem(alps::params & self, nb::object const & key_obj, nb::handle value) {
    std::string const key = nb::cast<std::string>(nb::str(key_obj));
    std::string const text = pyalps::param_text(value, key);
    if (self.exists(key)) {
        alps::params parsed;
        parsed.supply(key, text);
        if (self[key].apply_visitor(parse_as{{}, parsed, key})) {
            self[key] = parsed[key];
            return;
        }
    }
    if (self.exists(key) || self.supplied(key)) {
        alps::params copy = without(self, key);
        copy.supply(key, text);
        swap(self, copy);
    } else
        self.supply(key, text);
}
void params_delitem(alps::params & self, nb::object const & key_obj) {
    std::string const key = nb::cast<std::string>(nb::str(key_obj));
    if (!params_contains(self, key))
        throw nb::key_error(key.c_str());
    alps::params copy = without(self, key);
    swap(self, copy);
}
std::vector<std::string> params_keys(alps::params const & self) {
    std::set<std::string> keys;
    for (auto const & entry : self)
        keys.insert(entry.first);
    for (auto const & supplied : self.supplied_text())
        keys.insert(supplied.first);
    return std::vector<std::string>(keys.begin(), keys.end());
}
nb::object value_or_default(alps::params const & self, nb::object const & key, nb::object const & dflt) {
    return params_contains(self, nb::cast<std::string>(nb::str(key))) ? params_getitem(self, key) : dflt;
}
void params_load(alps::params & self, alps::hdf5::archive & ar, std::string const & path) {
    ar[path] >> self;
}
std::string params_print(alps::params const & self) {
    std::stringstream ss;
    ss << self;
    return ss.str();
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
                 new (self) alps::params(alps::make_parameters_from_file(filename));
             },
             nb::arg("filename"))
        .def("__init__", [](alps::params * self, alps::hdf5::archive & ar, std::string const & path) {
                 new (self) alps::params();
                 params_load(*self, ar, path);
             },
             nb::arg("archive"),
             nb::arg("path") = std::string("/parameters"))
        .def("__len__",      [](alps::params const & self) { return params_keys(self).size(); })
        .def("__deepcopy__", [](alps::params const & self, nb::handle) { return alps::params(self); })
        .def("__getitem__",  &params_getitem)
        .def("__setitem__",  &params_setitem, nb::arg("key"), nb::arg("value").none())
        .def("__delitem__",  &params_delitem)
        .def("__contains__", [](alps::params const & self, nb::object const & key) {
                                 return params_contains(self, nb::cast<std::string>(nb::str(key)));
                             })
        .def("__iter__",     [](alps::params const & self) {
                                 // A snapshot: assignments may rebuild the object.
                                 return nb::cast(params_keys(self)).attr("__iter__")();
                             })
        .def("__str__",      &params_print)
        .def("valueOrDefault", &value_or_default)
        .def("supplied",     &alps::params::supplied, nb::arg("key"),
             "True if the parameter was supplied, by a file, a dict or an assignment")
        .def("save",         &alps::params::save)
        .def("load",         &params_load,
             nb::arg("archive"),
             nb::arg("path") = std::string("/parameters"));
    pyalps::mark_archive_savable(m.attr("params"));
}
