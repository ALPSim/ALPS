// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// CSR matrix view and the deflated operator the Lanczos driver hands to
// IETL. Include this before <ietl/lanczos.h>: IETL calls ietl::mult by
// qualified name inside its templates, so the overload for
// DeflatedOperator must be declared at the template's definition.

#pragma once

#include <boost/numeric/ublas/vector.hpp>

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace alps::ed::detail {

inline double conjugate(double value) { return value; }
inline std::complex<double> conjugate(std::complex<double> value) { return std::conj(value); }
inline double real_part(double value) { return value; }
inline double real_part(std::complex<double> value) { return value.real(); }

// Borrowed CSR buffers (row = output, column = input); entries may repeat and
// columns may be unordered, both are accumulated.
template<class Scalar> struct CsrView {
    std::int64_t dim = 0;
    std::int64_t const* offsets = nullptr;
    std::int64_t const* indices = nullptr;
    Scalar const* values = nullptr;

    template<class X, class Y> void apply(X const& x, Y& y) const {
        for (std::int64_t row = 0; row < dim; ++row) {
            Scalar sum{};
            for (auto p = offsets[row]; p < offsets[row + 1]; ++p)
                sum += values[p] * x[static_cast<std::size_t>(indices[p])];
            y[static_cast<std::size_t>(row)] = sum;
        }
    }
};

// <u|v> over the first n entries of two vector-like ranges.
template<class Scalar, class U, class V>
Scalar inner(U const& u, V const& v, std::size_t n) {
    Scalar sum{};
    for (std::size_t i = 0; i < n; ++i) sum += conjugate(Scalar(u[i])) * Scalar(v[i]);
    return sum;
}

// v -= sum_u u <u|v> for the locked orthonormal vectors u.
template<class Scalar, class V>
void project_out(std::vector<std::vector<Scalar>> const& locked, V& v, std::size_t n) {
    for (auto const& u : locked) {
        Scalar const overlap = inner<Scalar>(u, v, n);
        for (std::size_t i = 0; i < n; ++i) v[i] -= overlap * u[i];
    }
}

// y = (I - P)(sign*A - shift)(I - P) x, with P the projector onto the locked
// span. Locked directions have eigenvalue 0; with shift above the spectrum
// of sign*A every other eigenvalue is negative, so the lowest Ritz pairs of
// this operator are the lowest unlocked eigenpairs of sign*A, shifted.
template<class Scalar> struct DeflatedOperator {
    using Vector = boost::numeric::ublas::vector<Scalar>;
    CsrView<Scalar> const* matrix = nullptr;
    double sign = 1.;
    double shift = 0.;
    std::vector<std::vector<Scalar>> const* locked = nullptr;
    mutable std::vector<Scalar> scratch;

    void apply(Vector const& x, Vector& y) const {
        auto const n = static_cast<std::size_t>(matrix->dim);
        scratch.assign(x.begin(), x.end());
        if (!locked->empty()) project_out(*locked, scratch, n);
        matrix->apply(scratch, y);
        for (std::size_t i = 0; i < n; ++i) y[i] = sign * y[i] - shift * scratch[i];
        if (!locked->empty()) project_out(*locked, y, n);
    }
};

} // namespace alps::ed::detail

namespace ietl {
template<class Scalar>
void mult(alps::ed::detail::DeflatedOperator<Scalar> const& op,
          boost::numeric::ublas::vector<Scalar> const& x,
          boost::numeric::ublas::vector<Scalar>& y) {
    op.apply(x, y);
}
} // namespace ietl
