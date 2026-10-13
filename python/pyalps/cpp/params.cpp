// Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>
//               Matthias Troyer <troyer@comp-phys.org>
//               2026       by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <alps/hdf5/archive.hpp>
#include "archive_savable.hpp"
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
#include "dict_to_params.hpp"
namespace nb = nanobind;
namespace {
// Walk the paramvalue variant and wrap each native alternative as a
// nb::object. Called from __getitem__.
struct paramvalue_to_py_visitor : boost::static_visitor<nb::object> {
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
                nb::cast(value), nb::arg("dtype") = alps::python::numpy_dtype<T>::name);
    }
};
nb::object paramvalue_to_py(alps::detail::paramvalue const & pv) {
    if (pv.source()) {
        auto * value = static_cast<PyObject *>(pv.source()->object("python"));
        if (value)
            return nb::borrow<nb::object>(value);
        std::vector<alps::detail::paramvalue> elements;
        if (pv.source()->native_elements(elements)) {
            nb::list result;
            for (auto const & element : elements)
                result.append(paramvalue_to_py(element));
            return result;
        }
        return paramvalue_to_py(pv.source()->native_value());
    }
    return boost::apply_visitor(
        paramvalue_to_py_visitor(),
        static_cast<alps::detail::paramvalue_base const &>(pv));
}
// Retain the value through the shared binding-owned provider.
void params_setitem(alps::params & self, nb::object const & key_obj, nb::handle value) {
    pyalps::set_param_value(self, nb::cast<std::string>(nb::str(key_obj)), value);
}
nb::object params_getitem(alps::params & self, nb::object const & key_obj) {
    std::string key = nb::cast<std::string>(nb::str(key_obj));
    alps::detail::paramvalue const * value = self.find(key);
    if (!value)
        return nb::none();
    nb::object result = paramvalue_to_py(*value);
    // Materialize native checkpoint values once. Retain the object so later
    // mutations survive both subsequent Python lookups and C++ conversions.
    if (!value->source() || !value->source()->object("python"))
        pyalps::set_param_value(self, key, result);
    return result;
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
    alps::hdf5::archive reader(ar);
    reader.set_context(ar.complete_path(path));
    pyalps::enable_python_param_reader(self);
    self.load(reader);
}
std::string params_print(alps::params & self) {
    std::stringstream ss;
    ss << self;
    return ss.str();
}
// deepcopy support, including shared objects in the caller's memo.
alps::params params_deepcopy(alps::params & self, nb::handle memo) {
    nb::dict values;
    for (auto const & entry : self) {
        nb::str key(entry.first.c_str());
        values[key] = params_getitem(self, key);
    }
    return pyalps::params_from_dict(nb::cast<nb::dict>(
        nb::module_::import_("copy").attr("deepcopy")(values, memo)));
}
}  // namespace
NB_MODULE(pyngsparams_c, m) {
    nb::class_<alps::params>(m, "params")
        .def("__init__", [](alps::params * self) {
            new (self) alps::params(pyalps::params_from_dict(nb::dict()));
        })
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
                 pyalps::enable_python_param_reader(*self);
             },
             nb::arg("filename"))
        .def("__init__", [](alps::params * self, alps::hdf5::archive & ar, std::string const & path) {
                 alps::params loaded;
                 params_load(loaded, ar, path);
                 new (self) alps::params(loaded);
             },
             nb::arg("archive"),
             nb::arg("path") = std::string("/parameters"))
        .def("__len__",      [](alps::params const & self) { return self.size(); })
        .def("__deepcopy__", &params_deepcopy)
        .def("__getitem__",  &params_getitem)
        .def("__setitem__",  &params_setitem, nb::arg("key"), nb::arg("value").none())
        .def("__delitem__",  &params_delitem)
        .def("__contains__", &params_contains)
        .def("__iter__",     [](alps::params & self) {
                                 // Snapshot keys: params' native iterator
                                 // holds a vector iterator invalidated by
                                 // insertion or deletion, even while the
                                 // params object itself remains alive.
                                 nb::list keys;
                                 for (auto const & entry : self)
                                     keys.append(nb::cast(entry.first));
                                 return keys.attr("__iter__")();
                             })
        .def("__str__",      &params_print)
        .def("valueOrDefault", &value_or_default)
        .def("save",         &alps::params::save)
        .def("load",         &params_load,
             nb::arg("archive"),
             nb::arg("path") = std::string("/parameters"));
    pyalps::mark_archive_savable(m.attr("params"));
}
