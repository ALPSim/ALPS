// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// ALPS's addition to ALPSCore's params (see the README): the text that values
// from HDF5 groups without ALPSCore's attributes and from Python are supplied
// as, so that define<T>() parses them as it parses values from INI files.
#pragma once

#include <alps/hdf5/archive.hpp>
#include <alps/params.hpp>

#include <string>
#include <type_traits>
#include <vector>

namespace alps {
    namespace params_ns {

        /// Text of a floating-point number that parses back to the same value
        /** An integral value is written without a fraction or exponent, so
            that `24.0` also parses as an integer; a non-integral one has as
            many digits as it needs, and parsing it as an integer fails. */
        std::string to_text(double value);

        /// `1` or `0`, which parse both as `bool` and as integers
        inline std::string to_text(bool value) { return value ? "1" : "0"; }

        /// Decimal text of an integer
        template <typename T>
        typename std::enable_if<std::is_integral<T>::value && !std::is_same<T, bool>::value,
                                std::string>::type
        to_text(T value) { return std::to_string(value); }

        inline std::string to_text(const std::string& value) { return value; }

        /// The elements' text separated by commas, as vector parameters are written
        template <typename T>
        std::string to_text(const std::vector<T>& values) {
            std::string text;
            for (typename std::vector<T>::size_type i = 0; i < values.size(); ++i) {
                if (i != 0) text += ',';
                text += to_text(static_cast<T>(values[i]));
            }
            return text;
        }

        /// What supplied text reads as when no definition gives its type
        /** An integer that fits in a long, another number as a double,
            numbers separated by commas as a vector of long if each is such
            an integer and as a vector of double otherwise, or else the text
            itself. params::save() writes parameters that are supplied but
            not defined as these values, and pyalps reads them so. */
        dict_value supplied_value(const std::string& name, const std::string& text);

        /// Supplies each child of the archive's current group as text
        /** Scalars and one-dimensional arrays of numbers, booleans and
            strings are supplied with `to_text()`, and so are groups of
            scalars named `0`, `1`, ..., as Python lists are stored.
            `params::load()` calls this for a group without ALPSCore's
            attributes, as ALPS's other writers leave it.
            @throws std::runtime_error naming the parameter, for complex,
            multidimensional or other values that params cannot hold */
        void supply_plain_group(alps::hdf5::archive& ar, params& p);

    }
}
