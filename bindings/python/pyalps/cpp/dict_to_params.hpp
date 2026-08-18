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
#include <string>
#include <vector>
namespace pyalps {
namespace nb = nanobind;
// Store one Python value under `key`. paramvalue's only integral
// alternative is a 32-bit int and libalps static_casts wider integer
// types down to it, so out-of-range Python ints are rejected loudly
// here rather than truncated silently. List probes use exact element
// types first (convert=false) so integer lists round-trip as ints;
// mixed numeric lists without bools widen to double.
inline void set_param_value(alps::params & p, std::string const & key, nb::handle value) {
    if (value.is_none())
        throw nb::type_error(("cannot store None for parameter '" + key
            + "': params has no null type; delete the key instead").c_str());
    if (nb::isinstance<nb::bool_>(value)) {
        p[key] = nb::cast<bool>(value);
    } else if (nb::isinstance<nb::int_>(value)) {
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
        try { p[key] = nb::cast<std::vector<int>>(value, false); return; }
        catch (nb::cast_error const &) {}
        try { p[key] = nb::cast<std::vector<double>>(value, false); return; }
        catch (nb::cast_error const &) {}
        try { p[key] = nb::cast<std::vector<std::complex<double>>>(value, false); return; }
        catch (nb::cast_error const &) {}
        try { p[key] = nb::cast<std::vector<std::string>>(value, false); return; }
        catch (nb::cast_error const &) {}
        // mixed numeric content (e.g. [1, 2.5]) widens to double —
        // but never bools, which would silently become 0.0/1.0
        nb::object seq = nb::borrow<nb::object>(value);
        std::size_t const n = nb::len(seq);
        bool has_bool = false;
        for (std::size_t i = 0; i < n && !has_bool; ++i) {
            nb::object item = seq[i];
            has_bool = nb::isinstance<nb::bool_>(item);
        }
        if (!has_bool) {
            try { p[key] = nb::cast<std::vector<double>>(value); return; }
            catch (nb::cast_error const &) {}
            try { p[key] = nb::cast<std::vector<std::complex<double>>>(value); return; }
            catch (nb::cast_error const &) {}
        }
        throw nb::type_error(("unsupported list for parameter '" + key
            + "' (expected homogeneous numbers or strings)").c_str());
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
