// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Shared dense/Lanczos solve of one symmetry block held as canonical CSR.
// Factored out of the block solvers in skilledwolf/alpsx (branch
// agent/v3-execution-and-api-audit at f8ee6e98b, files
// src/solvers/diag/src/runtime/{projected_basis,spectrum}.cpp), which held
// their blocks as Eigen sparse matrices and solved them with Spectra.
#pragma once

#include <alps/ed/eigensolver.hpp>
#include <alps/ed/lanczos.hpp>
#include <alps/ed/sparse_matrix.hpp>
#include <alps/ed/spectrum.hpp>
#include "eigenstate_policy.hpp"

#include <complex>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace alps::ed::detail {

struct SolvedBlock {
    SpectrumBlock spectrum;
    std::vector<std::complex<double>> vectors;  // column-major, when requested
    std::vector<double> residuals;              // Lanczos only, per column
};

template<class Scalar>
void solve_block_as(SparseMatrix<std::complex<double>> const& matrix, std::string const& which,
                    SolvePlan const& plan, bool want_vectors, SolvedBlock& result) {
    auto const n = matrix.rows;
    auto const k = plan.count;
    auto const scalar = [](std::complex<double> z) -> Scalar {
        if constexpr (std::is_same_v<Scalar, double>) return z.real();
        else return z;
    };
    if (!plan.iterative) {
        std::vector<Scalar> dense;
        if (n > dense.max_size() / n) throw std::length_error("projected ED: dense block is too large");
        dense.resize(n * n);
        for (std::size_t row = 0; row < n; ++row)
            for (auto p = matrix.offsets[row]; p < matrix.offsets[row + 1]; ++p)
                dense[static_cast<std::size_t>(matrix.indices[static_cast<std::size_t>(p)]) * n + row]
                    = scalar(matrix.values[static_cast<std::size_t>(p)]);
        auto value = dense_eigensystem(std::move(dense), n, want_vectors, k, which);
        result.spectrum.eigenvalues = std::move(value.eigenvalues);
        result.vectors.assign(value.eigenvectors_column_major.begin(), value.eigenvectors_column_major.end());
        result.spectrum.method = "dense";
        result.spectrum.nconv = static_cast<int>(k);
        return;
    }
    std::vector<Scalar> values;
    values.reserve(matrix.values.size());
    for (auto z : matrix.values) values.push_back(scalar(z));
    auto value = lanczos_csr(matrix.offsets, matrix.indices, values, static_cast<std::int64_t>(n),
                             static_cast<int>(k), which, plan.tolerance, plan.max_iterations, want_vectors);
    result.spectrum.eigenvalues = std::move(value.eigenvalues);
    result.vectors.assign(value.eigenvectors.begin(), value.eigenvectors.end());
    result.residuals = std::move(value.residual_norms);
    result.spectrum.method = "lanczos";
    result.spectrum.nconv = value.nconv;
    result.spectrum.tolerance = value.tolerance;
    result.spectrum.max_iterations = value.max_iterations;
    result.spectrum.iterations = value.iterations;
}

// Callers validate the public controls and resolve the policy for this block.
// Real symmetric blocks (every stored imaginary part exactly zero) are solved
// in real arithmetic; the returned vectors are complex either way.
inline SolvedBlock solve_block(SparseMatrix<std::complex<double>> const& matrix, std::string const& which,
                               SolvePlan const& plan, bool want_vectors) {
    SolvedBlock result;
    result.spectrum.block_dimension = matrix.rows;
    if (matrix.rows == 0) {
        result.spectrum.method = "dense";
        return result;
    }
    bool real = true;
    for (auto z : matrix.values) real = real && z.imag() == 0.;
    if (real) solve_block_as<double>(matrix, which, plan, want_vectors, result);
    else solve_block_as<std::complex<double>>(matrix, which, plan, want_vectors, result);
    return result;
}

} // namespace alps::ed::detail
