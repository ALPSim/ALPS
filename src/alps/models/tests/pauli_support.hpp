// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Shared helpers of the Pauli export tests, adapted from the reference
// helpers of skilledwolf/alpsx (branch agent/v3-execution-and-api-audit at
// f8ee6e98b: src/model/tests/models/pauli_export.cpp, tests/alps/
// test_spin_pauli.py and test_fermion_pauli.py): the dense matrix of a Pauli
// operator, canonical anticommutation relations applied directly to
// occupation integers (independent of any Jordan-Wigner string), angular
// momentum matrices independent of the kit's site factories, and the
// injection of the physical product basis into the qubit register. Dense
// algebra is the row-major DenseMatrix of dense_e0.hpp; no Eigen.

#ifndef ALPS_TESTS_MODELS_PAULI_SUPPORT_HPP
#define ALPS_TESTS_MODELS_PAULI_SUPPORT_HPP

#include "dense_e0.hpp"
#include <alps/models/pauli.hpp>

#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace alps_test {

using Complex = std::complex<double>;

// 2^n x 2^n matrix of a Pauli operator; qubit q is bit q of the index and
// letter n-1-q of a word.
inline DenseMatrix pauli_matrix(alps::models::PauliOperator const& p) {
    std::size_t const dim = std::size_t{1} << p.num_qubits;
    DenseMatrix result(dim);
    for (auto const& [word, coefficient] : p.terms)
        for (std::size_t column = 0; column < dim; ++column) {
            std::size_t row = column;
            Complex weight = coefficient;
            for (int q = 0; q < p.num_qubits; ++q) {
                char const ch = word[static_cast<std::size_t>(p.num_qubits - 1 - q)];
                bool const bit = (column >> q) & 1;
                if (ch == 'X' || ch == 'Y') row ^= std::size_t{1} << q;
                if (ch == 'Z') weight *= bit ? -1. : 1.;
                if (ch == 'Y') weight *= Complex(0., bit ? -1. : 1.);
            }
            result(row, column) += weight;
        }
    return result;
}

// Annihilator of `mode` on 2^n occupation integers under the canonical
// anticommutation relations: mode m occupies bit mode_bits[m], and the sign
// counts the occupied modes before m in mode order.
inline DenseMatrix annihilator(int num_modes, int mode, std::vector<int> const& mode_bits) {
    std::size_t const dim = std::size_t{1} << num_modes;
    DenseMatrix result(dim);
    int const bit = mode_bits[static_cast<std::size_t>(mode)];
    for (std::size_t state = 0; state < dim; ++state)
        if ((state >> bit) & 1) {
            int occupied = 0;
            for (int before = 0; before < mode; ++before)
                occupied += static_cast<int>((state >> mode_bits[static_cast<std::size_t>(before)]) & 1);
            result(state ^ (std::size_t{1} << bit), state) = occupied % 2 ? -1. : 1.;
        }
    return result;
}
inline DenseMatrix annihilator(int num_modes, int mode) {
    std::vector<int> bits(static_cast<std::size_t>(num_modes));
    std::iota(bits.begin(), bits.end(), 0);
    return annihilator(num_modes, mode, bits);
}

inline DenseMatrix multiply(DenseMatrix const& a, DenseMatrix const& b) {
    if (a.dim != b.dim) throw std::invalid_argument("multiply: dimension mismatch");
    DenseMatrix out(a.dim);
    for (std::size_t r = 0; r < a.dim; ++r)
        for (std::size_t k = 0; k < a.dim; ++k) {
            auto const x = a(r, k);
            if (x == 0.) continue;
            for (std::size_t c = 0; c < a.dim; ++c) out(r, c) += x * b(k, c);
        }
    return out;
}
inline DenseMatrix adjoint(DenseMatrix const& a) {
    DenseMatrix out(a.dim);
    for (std::size_t r = 0; r < a.dim; ++r)
        for (std::size_t c = 0; c < a.dim; ++c) out(c, r) = std::conj(a(r, c));
    return out;
}
inline DenseMatrix add(DenseMatrix a, DenseMatrix const& b, Complex scale = 1.) {
    if (a.dim != b.dim) throw std::invalid_argument("add: dimension mismatch");
    for (std::size_t i = 0; i < a.data.size(); ++i) a.data[i] += scale * b.data[i];
    return a;
}
inline DenseMatrix scale(DenseMatrix a, Complex factor) {
    for (auto& z : a.data) z *= factor;
    return a;
}
inline double max_abs(DenseMatrix const& a) {
    double m = 0.;
    for (auto z : a.data) m = std::max(m, std::abs(z));
    return m;
}
inline double max_abs_difference(DenseMatrix const& a, DenseMatrix const& b) { return max_abs(add(a, b, -1.)); }
inline std::vector<Complex> matvec(DenseMatrix const& a, std::vector<Complex> const& v) {
    if (v.size() != a.dim) throw std::invalid_argument("matvec: dimension mismatch");
    std::vector<Complex> out(a.dim);
    for (std::size_t r = 0; r < a.dim; ++r)
        for (std::size_t c = 0; c < a.dim; ++c) out[r] += a(r, c) * v[c];
    return out;
}
inline Complex vdot(std::vector<Complex> const& a, std::vector<Complex> const& b) {
    Complex out = 0.;
    for (std::size_t i = 0; i < a.size(); ++i) out += std::conj(a[i]) * b[i];
    return out;
}

// Spin-S matrices (x, y, z) with m = S, S-1, ..., -S, independent of the kit.
struct SpinMatrices { DenseMatrix x, y, z; };
inline SpinMatrices spin_matrices(double S) {
    std::size_t const d = static_cast<std::size_t>(std::lround(2. * S)) + 1;
    DenseMatrix raising(d), z(d);
    for (std::size_t column = 0; column < d; ++column) {
        double const m = S - static_cast<double>(column);
        z(column, column) = m;
        if (column > 0) raising(column - 1, column) = std::sqrt((S - m) * (S + m + 1.));
    }
    auto const lowering = adjoint(raising);
    return {scale(add(raising, lowering), .5), scale(add(raising, lowering, -1.), Complex(0., -.5)), z};
}

// Register indices of the physical product states in site-0-most-significant
// order (the order of a Kronecker product over the sites).
inline std::vector<std::uint64_t> physical_indices(alps::models::QubitEncoding const& encoding) {
    auto const& dims = encoding.site_dimensions();
    std::size_t total = 1;
    for (int d : dims) total *= static_cast<std::size_t>(d);
    std::vector<std::uint64_t> out;
    std::vector<int> indices(dims.size());
    for (std::size_t state = 0; state < total; ++state) {
        std::size_t rest = state;
        for (std::size_t s = dims.size(); s-- > 0;) {
            indices[s] = static_cast<int>(rest % static_cast<std::size_t>(dims[s]));
            rest /= static_cast<std::size_t>(dims[s]);
        }
        out.push_back(encoding.state_index(indices));
    }
    return out;
}

// The physical block A[indices, indices] of a register matrix.
inline DenseMatrix restrict_to(DenseMatrix const& a, std::vector<std::uint64_t> const& indices) {
    DenseMatrix out(indices.size());
    for (std::size_t r = 0; r < indices.size(); ++r)
        for (std::size_t c = 0; c < indices.size(); ++c)
            out(r, c) = a(static_cast<std::size_t>(indices[r]), static_cast<std::size_t>(indices[c]));
    return out;
}

// max |(A V - V B)| for the injection V of the physical basis: A acts on
// physical states exactly as B does (the columns of A V are physical states).
inline double action_distance(DenseMatrix const& a, std::vector<std::uint64_t> const& indices, DenseMatrix const& b) {
    if (b.dim != indices.size()) throw std::invalid_argument("action_distance: dimension mismatch");
    double m = 0.;
    for (std::size_t c = 0; c < indices.size(); ++c)
        for (std::size_t r = 0; r < a.dim; ++r) {
            Complex av = a(r, static_cast<std::size_t>(indices[c]));
            Complex vb = 0.;
            for (std::size_t k = 0; k < indices.size(); ++k)
                if (static_cast<std::size_t>(indices[k]) == r) vb = b(k, c);
            m = std::max(m, std::abs(av - vb));
        }
    return m;
}

// Largest column norm of A restricted to the given columns (zero when A
// annihilates every one of them).
inline double max_column_norm(DenseMatrix const& a, std::vector<std::uint64_t> const& columns) {
    double m = 0.;
    for (auto c : columns) {
        double sum = 0.;
        for (std::size_t r = 0; r < a.dim; ++r) sum += std::norm(a(r, static_cast<std::size_t>(c)));
        m = std::max(m, std::sqrt(sum));
    }
    return m;
}

} // namespace alps_test

#endif // ALPS_TESTS_MODELS_PAULI_SUPPORT_HPP
