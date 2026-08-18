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
inline bool is_bool_like(PyObject * raw) {
    // plain bool, or a numpy bool scalar (numpy.bool_ / numpy.bool),
    // which does NOT subclass bool and would otherwise slip through
    // the numeric ladder as 0.0/1.0
    return PyBool_Check(raw)
        || std::strncmp(Py_TYPE(raw)->tp_name, "numpy.bool", 10) == 0;
}
} // namespace detail
// Store one Python value under `key`. paramvalue's only integral
// alternative is a 32-bit int and libalps static_casts wider integer
// types down to it, so out-of-range integers are rejected loudly here
// — for scalars and inside lists alike — rather than truncated or
// silently widened to double. List probes use exact element types
// first (convert=false) so integer lists round-trip as ints; mixed
// numeric lists without bools widen to double.
inline void set_param_value(alps::params & p, std::string const & key, nb::handle value) {
    if (value.is_none())
        throw nb::type_error(("cannot store None for parameter '" + key
            + "': params has no null type; delete the key instead").c_str());
    if (detail::is_bool_like(value.ptr())) {
        // PyObject_IsTrue rather than nb::cast<bool>: the caster does
        // not convert numpy bool scalars
        int const truth = PyObject_IsTrue(value.ptr());
        if (truth < 0)
            throw nb::python_error();
        p[key] = (truth == 1);
    } else if (nb::isinstance<nb::int_>(value) || PyIndex_Check(value.ptr())) {
        // PyIndex_Check admits numpy integer scalars (np.int64 etc.),
        // which don't subclass int the way np.float64 subclasses float
        try {
            p[key] = nb::cast<int>(value);
        } catch (nb::cast_error const &) {
            throw nb::type_error(("parameter '" + key
                + "' does not fit params' 32-bit integer type").c_str());
        }
    } else if (nb::isinstance<nb::float_>(value)) {
        p[key] = nb::cast<double>(value);
    } else if (PyComplex_Check(value.ptr())) {
        p[key] = nb::cast<std::complex<double>>(value);
    } else if (nb::isinstance<nb::str>(value)) {
        p[key] = nb::cast<std::string>(value);
    } else if (nb::isinstance<nb::list>(value) || nb::isinstance<nb::tuple>(value)) {
        // One pre-scan enforcing the loud-failure policies explicitly,
        // independent of caster conversion behaviour: bools never
        // coerce to numbers, and oversized integers raise exactly like
        // the scalar arm instead of widening to double (which would
        // corrupt values beyond 2^53).
        std::size_t const length = nb::len(value);
        bool has_bool = false;
        for (std::size_t i = 0; i < length; ++i) {
            nb::object item = value[i];
            PyObject * raw = item.ptr();
            if (detail::is_bool_like(raw)) {
                has_bool = true;
            } else if (PyLong_Check(raw) || PyIndex_Check(raw)) {
                // PyNumber_Index covers numpy integer scalars too —
                // they are not PyLong subclasses but must obey the
                // same 32-bit range policy
                PyObject * as_long = PyNumber_Index(raw);
                if (!as_long) {
                    PyErr_Clear();
                    continue;
                }
                int overflow = 0;
                long long v = PyLong_AsLongLongAndOverflow(as_long, &overflow);
                Py_DECREF(as_long);
                if (overflow
                    || v < std::numeric_limits<int>::min()
                    || v > std::numeric_limits<int>::max())
                    throw nb::type_error(("parameter '" + key
                        + "' contains an integer that does not fit params'"
                          " 32-bit integer type").c_str());
            }
        }
        if (!has_bool) {
            try { p[key] = nb::cast<std::vector<int>>(value, false); return; }
            catch (nb::cast_error const &) {}
            try { p[key] = nb::cast<std::vector<double>>(value, false); return; }
            catch (nb::cast_error const &) {}
            try { p[key] = nb::cast<std::vector<std::complex<double>>>(value, false); return; }
            catch (nb::cast_error const &) {}
            try { p[key] = nb::cast<std::vector<std::string>>(value, false); return; }
            catch (nb::cast_error const &) {}
            // numpy integer scalars satisfy the convert=true int
            // caster via __index__ (floats don't), keeping
            // [np.int64(8)] consistent with the scalar np.int64 rung;
            // the pre-scan above already range-checked every element
            try { p[key] = nb::cast<std::vector<int>>(value); return; }
            catch (nb::cast_error const &) {}
            // mixed numeric content (e.g. [1, 2.5] or numpy floats)
            // widens to double / complex
            try { p[key] = nb::cast<std::vector<double>>(value); return; }
            catch (nb::cast_error const &) {}
            try { p[key] = nb::cast<std::vector<std::complex<double>>>(value); return; }
            catch (nb::cast_error const &) {}
        }
        throw nb::type_error(("unsupported list for parameter '" + key
            + "' (expected homogeneous numbers or strings; bools are not"
              " a parameter list type)").c_str());
    } else {
        throw nb::type_error(("unsupported type for parameter '" + key
            + "' (expected bool/int/float/complex/str or a list of those)").c_str());
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
