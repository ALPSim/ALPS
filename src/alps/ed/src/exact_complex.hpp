// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/src/runtime/exact_complex.hpp.
// Adaptation: C++17 (no defaulted comparisons), namespace alps::ed::detail.
#pragma once

#include <boost/multiprecision/cpp_int.hpp>
#include <complex>
#include <utility>

namespace alps::ed::detail {

// Finite binary floating-point inputs convert exactly. This private value is
// used to certify forbidden charge transitions, not to assemble or solve the
// allowed Hamiltonian blocks.
struct ExactComplex {
    using Rational = boost::multiprecision::cpp_rational;
    Rational real{}, imag{};

    ExactComplex() = default;
    explicit ExactComplex(std::complex<double> value) : real(value.real()), imag(value.imag()) {}

    bool zero() const { return real == 0 && imag == 0; }
    bool operator==(ExactComplex const& other) const { return real == other.real && imag == other.imag; }
    bool operator!=(ExactComplex const& other) const { return !(*this == other); }

    ExactComplex& operator+=(ExactComplex const& other) {
        real += other.real;
        imag += other.imag;
        return *this;
    }

    ExactComplex& operator*=(ExactComplex const& other) {
        Rational next_real = real * other.real - imag * other.imag;
        Rational next_imag = real * other.imag + imag * other.real;
        real = std::move(next_real);
        imag = std::move(next_imag);
        return *this;
    }
};

} // namespace alps::ed::detail
