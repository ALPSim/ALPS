// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// ALPS's addition to ALPSCore's params; see text_values.hpp and the README.

#include <alps/params/text_values.hpp>

#include <alps/hdf5/vector.hpp>

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>

namespace alps {
    namespace params_ns {

        std::string to_text(double value) {
            char buffer[64];
            if (std::isfinite(value) && value == std::floor(value) && std::fabs(value) < 1e18) {
                std::snprintf(buffer, sizeof buffer, "%.0f", value);
                return buffer;
            }
            for (int digits = 15; digits <= 17; ++digits) {
                std::snprintf(buffer, sizeof buffer, "%.*g", digits, value);
                if (!std::isfinite(value) || std::strtod(buffer, nullptr) == value) break;
            }
            return buffer;
        }

        namespace {

            // 1 for an integer that fits in a long, 2 for another number, 0 for
            // text that is not a number as a whole (hexadecimal is not one)
            int parse_number(const std::string& text, long& integer, double& real) {
                if (text.empty() || text.find_first_of("xX") != std::string::npos) return 0;
                const char* begin = text.c_str();
                char* end = nullptr;
                errno = 0;
                integer = std::strtol(begin, &end, 10);
                if (*end == '\0' && errno == 0) {
                    real = static_cast<double>(integer);
                    return 1;
                }
                real = std::strtod(begin, &end);
                return *end == '\0' ? 2 : 0;
            }

        } // namespace

        dict_value supplied_value(const std::string& name, const std::string& text) {
            dict_value value(name);
            long integer = 0;
            double real = 0.0;
            if (text.find(',') == std::string::npos) {
                switch (parse_number(text, integer, real)) {
                case 1: value = integer; break;
                case 2: value = real; break;
                default: value = text;
                }
                return value;
            }
            std::vector<long> integers;
            std::vector<double> reals;
            bool all_integers = true;
            for (std::string::size_type start = 0;;) {
                const std::string::size_type comma = text.find(',', start);
                const int kind = parse_number(
                    text.substr(start, comma == std::string::npos ? comma : comma - start), integer, real);
                if (kind == 0) {
                    value = text;
                    return value;
                }
                all_integers = all_integers && kind == 1;
                integers.push_back(integer);
                reals.push_back(real);
                if (comma == std::string::npos) break;
                start = comma + 1;
            }
            if (all_integers) value = integers;
            else value = reals;
            return value;
        }

        namespace {

            std::string error_prefix(const std::string& name) {
                return "parameter '" + name + "' in an HDF5 group";
            }

            // Reads `path` as a T or a vector of T if the archive stores that type
            template <typename T>
            bool read_as(alps::hdf5::archive& ar, const std::string& path, bool scalar, std::string& text) {
                if (!ar.is_datatype<T>(path)) return false;
                if (scalar) {
                    T value;
                    ar[path] >> value;
                    text = to_text(value);
                } else {
                    std::vector<T> values;
                    ar[path] >> values;
                    text = to_text(values);
                }
                return true;
            }

            // Floating-point types are read as double, signed char (also the
            // archive's bool) as int, and the other integers at their width.
            template <typename T>
            bool read_widened(alps::hdf5::archive& ar, const std::string& path, bool scalar, std::string& text) {
                typedef typename std::conditional<std::is_floating_point<T>::value, double,
                    typename std::conditional<std::is_same<T, signed char>::value, int, T>::type>::type read_type;
                if (!ar.is_datatype<T>(path)) return false;
                if (scalar) {
                    read_type value;
                    ar[path] >> value;
                    text = to_text(value);
                } else {
                    std::vector<read_type> values;
                    ar[path] >> values;
                    text = to_text(values);
                }
                return true;
            }

            std::string data_text(alps::hdf5::archive& ar, const std::string& path, const std::string& name) {
                if (ar.is_complex(path)) {
                    throw std::runtime_error(error_prefix(name) + " is complex, which params cannot hold");
                }
                const bool scalar = ar.is_scalar(path);
                if (!scalar && ar.dimensions(path) != 1) {
                    throw std::runtime_error(error_prefix(name) + " is a multidimensional array,"
                                             " which params cannot hold");
                }
                std::string text;
                if (read_as<std::string>(ar, path, scalar, text)
                    || read_widened<double>(ar, path, scalar, text)
                    || read_widened<float>(ar, path, scalar, text)
                    || read_widened<long double>(ar, path, scalar, text)
                    || read_widened<signed char>(ar, path, scalar, text)
                    || read_as<unsigned char>(ar, path, scalar, text)
                    || read_as<short>(ar, path, scalar, text)
                    || read_as<unsigned short>(ar, path, scalar, text)
                    || read_as<int>(ar, path, scalar, text)
                    || read_as<unsigned int>(ar, path, scalar, text)
                    || read_as<long>(ar, path, scalar, text)
                    || read_as<unsigned long>(ar, path, scalar, text)
                    || read_as<long long>(ar, path, scalar, text)
                    || read_as<unsigned long long>(ar, path, scalar, text)) {
                    return text;
                }
                throw std::runtime_error(error_prefix(name) + " has a datatype params cannot hold");
            }

        } // namespace

        void supply_plain_group(alps::hdf5::archive& ar, params& p) {
            for (const std::string& name : ar.list_children(ar.get_context())) {
                if (!ar.is_group(name)) {
                    p.supply(name, data_text(ar, name, name));
                    continue;
                }
                // A Python list of mixed or Boolean elements: scalars 0, 1, ...
                const std::vector<std::string> elements = ar.list_children(name);
                std::string text;
                for (std::size_t i = 0; i < elements.size(); ++i) {
                    const std::string element = name + "/" + std::to_string(i);
                    if (!ar.is_data(element) || !ar.is_scalar(element)) {
                        throw std::runtime_error(error_prefix(name) + " is a group that is not a list of scalars");
                    }
                    if (i != 0) text += ',';
                    text += data_text(ar, element, name);
                }
                p.supply(name, text);
            }
        }

    }
}
