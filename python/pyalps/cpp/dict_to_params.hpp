// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The single Python→alps::params conversion, shared by pyngsparams_c
// (__setitem__ / dict ctor), pyngsbase_c (mcbase ctor) and the application
// modules (maxent_c / cthyb / ctint), so every module ingests parameters
// identically. Python values are supplied as text, as an INI file or HDF5
// group would supply them (alps::params::supply), and each program's
// define<T>() parses the text with its own type; see
// src/alps/params/README.md.
#ifndef PYALPS_DICT_TO_PARAMS_HPP
#define PYALPS_DICT_TO_PARAMS_HPP
#include "numpy_compat.hpp"
#include <alps/params.hpp>
#include <alps/params/text_values.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <cstring>
#include <limits>
#include <string>
#include <vector>
namespace pyalps {
namespace nb = nanobind;
namespace detail {
enum class scalar_kind {
    unsupported,
    boolean,
    integer,
    real,
    complex,
    string
};

inline bool is_bool_like(PyObject * raw) {
    // plain bool, or a numpy bool scalar (numpy.bool_ / numpy.bool),
    // which does NOT subclass bool and would otherwise slip through
    // the numeric ladder as 0.0/1.0
    return PyBool_Check(raw)
        || alps::python::type_name(raw).compare(0, 10, "numpy.bool") == 0;
}
inline bool is_numpy_array(nb::handle value) {
    // isinstance, rather than an exact tp_name comparison, keeps ndarray
    // subclasses (for example an unmasked numpy.ma.MaskedArray) on the same
    // native-copy path. Reuse the process-lifetime module handle shared by
    // the other NumPy conversion helpers.
    return nb::isinstance(
        value, alps::python::numpy_module().attr("ndarray"));
}

inline char numpy_scalar_kind(nb::handle value) {
    nb::handle numpy = alps::python::numpy_module();
    if (!nb::isinstance(value, numpy.attr("generic")))
        return '\0';
    std::string const kind = nb::cast<std::string>(
        value.attr("dtype").attr("kind"));
    return kind.empty() ? '\0' : kind.front();
}

inline scalar_kind classify_scalar(nb::handle value) {
    if (is_bool_like(value.ptr()))
        return scalar_kind::boolean;

    // NumPy extended scalars (float16/32/longdouble and
    // complex64/clongdouble) are not consistently Python float/complex
    // subclasses. Inspect dtype.kind before consulting Python's protocols so
    // a complex scalar can never pass through __float__ and lose its
    // imaginary component.
    switch (numpy_scalar_kind(value)) {
        case 'b': return scalar_kind::boolean;
        case 'i':
        case 'u': return scalar_kind::integer;
        case 'f': return scalar_kind::real;
        case 'c': return scalar_kind::complex;
        case 'S':
        case 'U': return scalar_kind::string;
        case '\0': break;
        default: return scalar_kind::unsupported;
    }

    if (PyLong_Check(value.ptr()) || PyIndex_Check(value.ptr()))
        return scalar_kind::integer;
    if (PyFloat_Check(value.ptr()))
        return scalar_kind::real;
    if (PyComplex_Check(value.ptr()))
        return scalar_kind::complex;
    if (PyUnicode_Check(value.ptr()) || PyBytes_Check(value.ptr()))
        return scalar_kind::string;
    return scalar_kind::unsupported;
}

inline std::string string_value(nb::handle value) {
    if (PyUnicode_Check(value.ptr()))
        return nb::cast<std::string>(value);

    // Python/NumPy byte strings map to ALPS' native UTF-8 std::string. This
    // accepts the common fixed-width NumPy "S" dtype without creating an
    // opaque-object escape hatch; invalid UTF-8 remains a loud error because
    // nanobind must also decode the value when returning it to Python.
    char * bytes = nullptr;
    Py_ssize_t size = 0;
    if (PyBytes_AsStringAndSize(value.ptr(), &bytes, &size) != 0)
        throw nb::python_error();
    nb::object decoded = nb::steal<nb::object>(
        PyUnicode_DecodeUTF8(bytes, size, "strict"));
    if (!decoded.is_valid())
        throw nb::python_error();
    return nb::cast<std::string>(decoded);
}
} // namespace detail
// The text of a Python value: alps::params_ns::to_text's rules, with exact
// decimal integers of any size (define<T>() range-checks them) and
// comma-separated sequences. Values params cannot hold are a TypeError naming
// the parameter.
inline std::string param_text(nb::handle value, std::string const & key) {
    auto unsupported = [&key](std::string const & what) {
        return nb::type_error(("parameter '" + key + "' " + what).c_str());
    };
    if (value.is_none())
        throw unsupported("is None; params has no null type, delete the key instead");
    if (detail::is_numpy_array(value)) {
        std::size_t const ndim = nb::cast<std::size_t>(value.attr("ndim"));
        if (ndim == 0)
            return param_text(value.attr("item")(), key);
        if (ndim != 1)
            throw unsupported("is a multidimensional numpy array; params holds only scalars"
                              " and one-dimensional sequences");
        return param_text(value.attr("tolist")(), key);
    }
    switch (detail::classify_scalar(value)) {
        case detail::scalar_kind::boolean: {
            // PyObject_IsTrue rather than nb::cast<bool>: the caster does
            // not convert numpy bool scalars
            int const truth = PyObject_IsTrue(value.ptr());
            if (truth < 0)
                throw nb::python_error();
            return alps::params_ns::to_text(truth == 1);
        }
        case detail::scalar_kind::integer: {
            nb::object integer = nb::steal<nb::object>(PyNumber_Index(value.ptr()));
            if (!integer.is_valid())
                throw nb::python_error();
            return nb::cast<std::string>(nb::str(integer));
        }
        case detail::scalar_kind::real: {
            double const number = PyFloat_AsDouble(value.ptr());
            if (PyErr_Occurred())
                throw nb::python_error();
            return alps::params_ns::to_text(number);
        }
        case detail::scalar_kind::complex:
            throw unsupported("is complex, which params cannot hold");
        case detail::scalar_kind::string:
            return detail::string_value(value);
        case detail::scalar_kind::unsupported:
            break;
    }
    if (nb::isinstance<nb::list>(value) || nb::isinstance<nb::tuple>(value)) {
        std::string text;
        for (std::size_t i = 0; i < nb::len(value); ++i) {
            nb::object item = value[i];
            if (nb::isinstance<nb::list>(item) || nb::isinstance<nb::tuple>(item)
                || detail::is_numpy_array(item))
                throw unsupported("is a nested sequence; params holds only one-dimensional sequences");
            std::string const element = param_text(item, key);
            if (element.find(',') != std::string::npos)
                throw unsupported("is a sequence with an element containing a comma,"
                                  " which separates the elements");
            if (i) text += ',';
            text += element;
        }
        return text;
    }
    throw unsupported("has an unsupported type (expected bool/int/float/str, a one-dimensional"
                      " numpy array, or a sequence of those scalar types)");
}

inline alps::params params_from_dict(nb::dict const & values) {
    alps::params result;
    for (auto item : values) {
        std::string const key = nb::cast<std::string>(nb::str(item.first));
        result.supply(key, param_text(item.second, key));
    }
    return result;
}
} // namespace pyalps
#endif
