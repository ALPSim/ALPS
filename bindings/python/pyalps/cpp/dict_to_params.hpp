// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
#ifndef PYALPS_DICT_TO_PARAMS_HPP
#define PYALPS_DICT_TO_PARAMS_HPP
#include <alps/ngs/params.hpp>
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>
#include <string>
#include <vector>
namespace pyalps {
namespace nb = nanobind;
inline alps::params params_from_dict(nb::dict const & values) {
    alps::params result;
    for (auto item : values) {
        std::string key = nb::cast<std::string>(nb::str(item.first));
        nb::handle value = item.second;
        if (nb::isinstance<nb::bool_>(value))
            result[key] = nb::cast<bool>(value);
        else if (nb::isinstance<nb::int_>(value))
            result[key] = nb::cast<long>(value);
        else if (nb::isinstance<nb::float_>(value))
            result[key] = nb::cast<double>(value);
        else if (nb::isinstance<nb::str>(value))
            result[key] = nb::cast<std::string>(value);
        else if (nb::isinstance<nb::list>(value) || nb::isinstance<nb::tuple>(value))
            result[key] = nb::cast<std::vector<double>>(value);
        else
            throw nb::type_error(("unsupported parameter type for '" + key + "'").c_str());
    }
    return result;
}
} // namespace pyalps
#endif
