// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The single Python→alps::params conversion ladder, shared by
// pyngsparams_c (__setitem__ / dict ctor), pyngsbase_c (mcbase ctor)
// and the application modules (maxent_c / cthyb / ctint), so every
// module ingests parameters identically.
#ifndef PYALPS_DICT_TO_PARAMS_HPP
#define PYALPS_DICT_TO_PARAMS_HPP
#include "numpy_compat.hpp"
#include <alps/ngs/params.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <complex>
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
        || std::strncmp(Py_TYPE(raw)->tp_name, "numpy.bool", 10) == 0;
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

inline int integer_value(nb::handle value, std::string const & key) {
    PyObject * indexed = PyNumber_Index(value.ptr());
    if (!indexed)
        throw nb::python_error();
    int overflow = 0;
    long long const converted = PyLong_AsLongLongAndOverflow(indexed, &overflow);
    Py_DECREF(indexed);
    if (PyErr_Occurred()) {
        PyErr_Clear();
        overflow = 1;
    }
    if (overflow
        || converted < std::numeric_limits<int>::min()
        || converted > std::numeric_limits<int>::max())
        throw nb::type_error(("parameter '" + key
            + "' contains an integer that does not fit params'"
              " 32-bit integer type").c_str());
    return static_cast<int>(converted);
}

inline double real_value(nb::handle value, std::string const & key) {
    if (classify_scalar(value) == scalar_kind::integer)
        return static_cast<double>(integer_value(value, key));
    double const converted = PyFloat_AsDouble(value.ptr());
    if (PyErr_Occurred())
        throw nb::python_error();
    return converted;
}

inline std::complex<double> complex_value(nb::handle value,
                                          std::string const & key) {
    scalar_kind const kind = classify_scalar(value);
    if (kind == scalar_kind::integer || kind == scalar_kind::real)
        return std::complex<double>(real_value(value, key), 0.0);
    Py_complex const converted = PyComplex_AsCComplex(value.ptr());
    if (PyErr_Occurred())
        throw nb::python_error();
    return std::complex<double>(converted.real, converted.imag);
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
// Store one Python value under `key`. paramvalue's only integral
// alternative is a 32-bit int and libalps static_casts wider integer
// types down to it, so out-of-range integers are rejected loudly here
// — for scalars and inside lists alike — rather than truncated or
// silently widened to double. Sequence elements are classified before
// conversion so integers round-trip as ints, mixed real numerics widen to
// double, and complex values can never be coerced through a real-number path.
inline void set_param_value(alps::params & p, std::string const & key, nb::handle value) {
    if (value.is_none())
        throw nb::type_error(("cannot store None for parameter '" + key
            + "': params has no null type; delete the key instead").c_str());
    detail::scalar_kind const scalar_type = detail::classify_scalar(value);
    if (scalar_type == detail::scalar_kind::boolean) {
        // PyObject_IsTrue rather than nb::cast<bool>: the caster does
        // not convert numpy bool scalars
        int const truth = PyObject_IsTrue(value.ptr());
        if (truth < 0)
            throw nb::python_error();
        p[key] = (truth == 1);
    } else if (detail::is_numpy_array(value)) {
        // params is deliberately Python-object-free.  Convert the NumPy
        // value once at the boundary and store one of paramvalue's native
        // scalar/vector alternatives.  ALPS parameters are one-dimensional;
        // preserving an arbitrary N-D ndarray would require an object escape
        // hatch or a new tensor type in the C++ API.
        std::size_t const ndim = nb::cast<std::size_t>(value.attr("ndim"));
        if (ndim == 0) {
            set_param_value(p, key, value.attr("item")());
            return;
        }
        if (ndim != 1)
            throw nb::type_error(("parameter '" + key
                + "' is a multidimensional numpy array; params supports only scalars"
                  " and one-dimensional sequences").c_str());
        nb::object items = value.attr("tolist")();
        if (nb::len(items) == 0) {
            // With no elements the sequence ladder cannot infer a native
            // alternative. Preserve the ndarray's scalar family explicitly
            // so an empty bool array does not silently become vector<int>.
            std::string const kind = nb::cast<std::string>(
                value.attr("dtype").attr("kind"));
            if (kind == "b") p[key] = std::vector<bool>();
            else if (kind == "i" || kind == "u") p[key] = std::vector<int>();
            else if (kind == "f") p[key] = std::vector<double>();
            else if (kind == "c") p[key] = std::vector<std::complex<double>>();
            else if (kind == "S" || kind == "U") p[key] = std::vector<std::string>();
            else
                throw nb::type_error(("parameter '" + key
                    + "' is an empty numpy array with unsupported dtype kind '"
                    + kind + "'").c_str());
            return;
        }
        set_param_value(p, key, items);
        return;
    } else if (scalar_type == detail::scalar_kind::integer) {
        p[key] = detail::integer_value(value, key);
    } else if (scalar_type == detail::scalar_kind::real) {
        p[key] = detail::real_value(value, key);
    } else if (scalar_type == detail::scalar_kind::complex) {
        p[key] = detail::complex_value(value, key);
    } else if (scalar_type == detail::scalar_kind::string) {
        p[key] = detail::string_value(value);
    } else if (nb::isinstance<nb::list>(value) || nb::isinstance<nb::tuple>(value)) {
        // Classify before conversion. In particular, NumPy complex scalars
        // implement a warning-emitting __float__; blindly probing a
        // vector<double> caster first would discard their imaginary parts.
        std::size_t const length = nb::len(value);
        bool has_bool = false;
        bool has_integer = false;
        bool has_real = false;
        bool has_complex = false;
        bool has_string = false;
        for (std::size_t i = 0; i < length; ++i) {
            nb::object item = value[i];
            switch (detail::classify_scalar(item)) {
                case detail::scalar_kind::boolean: has_bool = true; break;
                case detail::scalar_kind::integer:
                    detail::integer_value(item, key); // range check now
                    has_integer = true;
                    break;
                case detail::scalar_kind::real: has_real = true; break;
                case detail::scalar_kind::complex: has_complex = true; break;
                case detail::scalar_kind::string: has_string = true; break;
                case detail::scalar_kind::unsupported:
                    throw nb::type_error(("unsupported element in parameter '"
                        + key + "' sequence").c_str());
            }
        }

        bool const has_non_bool = has_integer || has_real || has_complex || has_string;
        if (has_bool && !has_non_bool) {
            std::vector<bool> flags;
            flags.reserve(length);
            for (std::size_t i = 0; i < length; ++i) {
                int const truth = PyObject_IsTrue(value[i].ptr());
                if (truth < 0)
                    throw nb::python_error();
                flags.push_back(truth == 1);
            }
            p[key] = flags;
            return;
        }
        if (has_bool)
            throw nb::type_error(("unsupported sequence for parameter '" + key
                + "' (bools cannot be mixed with other element types)").c_str());
        if (has_string && (has_integer || has_real || has_complex))
            throw nb::type_error(("unsupported sequence for parameter '" + key
                + "' (strings cannot be mixed with numeric element types)").c_str());

        if (has_string) {
            std::vector<std::string> strings;
            strings.reserve(length);
            for (std::size_t i = 0; i < length; ++i)
                strings.push_back(detail::string_value(value[i]));
            p[key] = strings;
        } else if (has_complex) {
            std::vector<std::complex<double>> numbers;
            numbers.reserve(length);
            for (std::size_t i = 0; i < length; ++i)
                numbers.push_back(detail::complex_value(value[i], key));
            p[key] = numbers;
        } else if (has_real) {
            std::vector<double> numbers;
            numbers.reserve(length);
            for (std::size_t i = 0; i < length; ++i)
                numbers.push_back(detail::real_value(value[i], key));
            p[key] = numbers;
        } else {
            // An empty untyped Python sequence follows the historic native
            // ladder's first vector alternative (vector<int>).
            std::vector<int> numbers;
            numbers.reserve(length);
            for (std::size_t i = 0; i < length; ++i)
                numbers.push_back(detail::integer_value(value[i], key));
            p[key] = numbers;
        }
        return;
    } else {
        throw nb::type_error(("unsupported type for parameter '" + key
            + "' (expected bool/int/float/complex/str, a one-dimensional"
              " numpy array, or a sequence of those scalar types)").c_str());
    }
}
inline alps::params params_from_dict(nb::dict const & values) {
    alps::params result;
    for (auto item : values)
        set_param_value(result,
                        nb::cast<std::string>(nb::str(item.first)),
                        item.second);
    return result;
}
} // namespace pyalps
#endif
