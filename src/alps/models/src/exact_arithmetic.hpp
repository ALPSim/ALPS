// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/src/operators/exact_arithmetic.hpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.
#pragma once

#include <alps/models/site.hpp>
#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace alps::models::detail {

using Exact = boost::multiprecision::cpp_rational;

// Finite double inputs are dyadic rationals. Keep local products and sums
// exact until their explicit export boundary; this is not a public scalar type.
struct ExactComplex {
    Exact real{}, imag{};
    ExactComplex() = default;
    explicit ExactComplex(Complex value) : real(value.real()), imag(value.imag()) {}
    ExactComplex(Exact r, Exact i) : real(std::move(r)), imag(std::move(i)) {}
    bool zero() const { return real == 0 && imag == 0; }
    bool operator==(ExactComplex const& o) const { return real == o.real && imag == o.imag; }
    bool operator!=(ExactComplex const& o) const { return !(*this == o); }
    Exact norm_squared() const { return real * real + imag * imag; }
    ExactComplex& operator+=(ExactComplex const& other) {
        real += other.real; imag += other.imag;
        return *this;
    }
    ExactComplex operator-(ExactComplex const& other) const {
        return {real - other.real, imag - other.imag};
    }
    ExactComplex operator*(ExactComplex const& other) const {
        return {real * other.real - imag * other.imag,
                real * other.imag + imag * other.real};
    }
};

inline Exact scale_power_two(Exact value, std::int64_t exponent) {
    if (!exponent || value == 0) return value;
    if (exponent < -std::int64_t(std::numeric_limits<int>::max())
        || exponent > std::numeric_limits<int>::max())
        throw std::overflow_error("local product scale exceeds the supported exponent range");
    boost::multiprecision::cpp_int scale = 1;
    scale <<= static_cast<unsigned>(exponent < 0 ? -exponent : exponent);
    if (exponent > 0) value *= scale;
    else value /= scale;
    return value;
}

inline ExactComplex scale_power_two(ExactComplex const& value, std::int64_t exponent) {
    return {scale_power_two(value.real, exponent), scale_power_two(value.imag, exponent)};
}

inline bool export_component(Exact const& value, double& result) {
    result = value.convert_to<double>();
    if (!std::isfinite(result) || (result == 0. && value != 0)) return false;
    // An amplified partially rounded subnormal can have order-one relative
    // error. Move scale instead. Ordinary normal rounding remains explicit.
    return std::fpclassify(result) != FP_SUBNORMAL || Exact(result) == value;
}

inline bool export_complex(ExactComplex const& value, Complex& result) {
    double real{}, imag{};
    if (!export_component(value.real, real) || !export_component(value.imag, imag)) return false;
    result = {real, imag};
    return true;
}

struct ScaleInterval {
    std::int64_t lower = std::numeric_limits<std::int64_t>::lowest();
    std::int64_t upper = std::numeric_limits<std::int64_t>::max();

    // Exponents e for exporting value * 2^-e: finite normal rounding or an
    // exactly representable subnormal, with every nonzero component retained.
    void include(Exact const& value) {
        if (value == 0) return;
        boost::multiprecision::cpp_int n = numerator(value);
        if (n < 0) n = -n;
        auto const& d = denominator(value);
        auto const denominator_exponent = std::int64_t(boost::multiprecision::msb(d));
        auto const high = std::int64_t(boost::multiprecision::msb(n)) - denominator_exponent;
        auto const low = std::int64_t(boost::multiprecision::lsb(n)) - denominator_exponent;
        auto first = high - 1023;
        double converted{};
        if (!export_component(scale_power_two(value, -first), converted)) ++first;
        auto const last = std::max(high + 1022, low + 1074);
        lower = std::max(lower, first);
        upper = std::min(upper, last);
    }
    void include(ExactComplex const& value) { include(value.real); include(value.imag); }
};

inline std::int64_t checked_exponent_sum(std::int64_t a, std::int64_t b) {
    if ((b > 0 && a > std::numeric_limits<std::int64_t>::max() - b)
        || (b < 0 && a < std::numeric_limits<std::int64_t>::lowest() - b))
        throw std::overflow_error("local product scale sum exceeds the supported range");
    return a + b;
}

} // namespace alps::models::detail
