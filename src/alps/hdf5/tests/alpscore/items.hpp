// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// The values stored in archive_7146b9e1.h5. writer.cpp, built against
// ALPSCore, writes them; hdf5_alpscore_fixture.cpp, built against ALPS,
// reads them back. Both sides include only the standard library here.
#pragma once

#include <complex>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <valarray>
#include <vector>

#define ALPS_FIXTURE_ARG(...) __VA_ARGS__

// X(name, (type), value): stored at /data/<name>.
#define ALPS_FIXTURE_ITEMS(X) \
    X(int_scalar, (int), 42) \
    X(int64_scalar, (std::int64_t), std::int64_t(1) << 40) \
    X(uint_scalar, (unsigned), 7u) \
    X(int8_scalar, (signed char), static_cast<signed char>(-5)) \
    X(double_scalar, (double), 0.1) \
    X(float_scalar, (float), 0.25f) \
    X(bool_true, (bool), true) \
    X(bool_false, (bool), false) \
    X(string_scalar, (std::string), std::string("hello world")) \
    X(empty_string, (std::string), std::string()) \
    X(complex_scalar, (std::complex<double>), std::complex<double>(1.5, -2.5)) \
    X(vec_double, (std::vector<double>), (std::vector<double>{1.0, 2.0, 3.5})) \
    X(vec_int, (std::vector<int>), (std::vector<int>{1, -2, 3})) \
    X(vec_bool, (std::vector<bool>), (std::vector<bool>{true, false, true})) \
    X(vec_string, (std::vector<std::string>), (std::vector<std::string>{"a", "bc", ""})) \
    X(vec_complex, (std::vector<std::complex<double>>), (std::vector<std::complex<double>>{{1, 2}, {-3, 4}})) \
    X(vec_empty, (std::vector<double>), (std::vector<double>{})) \
    X(vec_vec_rect, (std::vector<std::vector<double>>), (std::vector<std::vector<double>>{{1, 2}, {3, 4}})) \
    X(vec_vec_ragged, (std::vector<std::vector<double>>), (std::vector<std::vector<double>>{{1}, {2, 3, 4}, {}})) \
    X(pair_di, (std::pair<double, int>), (std::pair<double, int>(2.5, 3))) \
    X(map_sd, (std::map<std::string, double>), (std::map<std::string, double>{{"x", 1.0}, {"y/z", 2.0}})) \
    X(map_empty, (std::map<std::string, double>), (std::map<std::string, double>{})) \
    X(valarray_d, (std::valarray<double>), (std::valarray<double>{1.0, 4.0, 9.0}))

namespace alps_fixture {
    template <class T> bool same(T const & a, T const & b) { return a == b; }
    inline bool same(std::valarray<double> const & a, std::valarray<double> const & b) {
        if (a.size() != b.size()) return false;
        for (std::size_t i = 0; i != a.size(); ++i)
            if (a[i] != b[i]) return false;
        return true;
    }
    // /data/multi_array is a 2x3 boost::multi_array<double, 2> holding i * 1.5.
    inline double multi_array_element(std::size_t i) { return static_cast<double>(i) * 1.5; }
}
