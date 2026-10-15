// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/include/alps/diag/eigensolver.hpp.
// Adaptation: C++17 (no std::span), <alps/ed/...> include names, ALPS_ED_DECL
// exports; the dense solver is LAPACK only (alpsx also used Eigen).
#pragma once

#include <alps/ed_export.h>

#include <complex>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace alps::ed {

// Model eigenstate controls shared by product problems. Automatic uses
// Lanczos above dimension 1000. A requested Lanczos solve uses dense
// diagonalization when the dimension does not exceed 8*count + 200, since
// the Lanczos driver needs Krylov room. Omitted Lanczos controls are the
// driver defaults (tolerance 1e-10, min(10*dimension, 1000) steps per run).
enum class EigenstateMethod { automatic, dense, lanczos };
struct EigenstateOptions {
    int count = 6;
    EigenstateMethod method = EigenstateMethod::automatic;
    std::string which = "SA";
    std::optional<double> tolerance;
    std::optional<int> max_iterations;
};

// Validate the shared public selection/convergence contract, including when
// a small block uses dense diagonalization. Omitted controls retain the
// calling workflow's documented defaults.
ALPS_ED_DECL void validate_eigensolver_controls(std::string const& which,
    std::optional<double> tolerance = {}, std::optional<int> max_iterations = {});

// Indices of the requested SA/LA/SM/LM subset, returned in ascending energy
// order. Magnitude ties use algebraic order; no particular basis is promised
// inside a degenerate eigenspace.
ALPS_ED_DECL std::vector<std::size_t> select_eigenvalues(std::vector<double> const& values,
    std::size_t count, std::string const& which);


template<class Scalar>
struct DenseEigenResult {
  std::vector<double> eigenvalues;
  std::vector<Scalar> eigenvectors_column_major;
};

// Consume an owned, finite column-major square matrix whose lower
// triangle defines the Hermitian operator. Outputs own their memory.
// Eigenvalues are ascending and eigenvectors are columns when requested.
// Omitted count returns the full spectrum; a positive count selects up to that
// many states using which, including when the selected states are not lowest.
ALPS_ED_DECL DenseEigenResult<double> dense_eigensystem(
    std::vector<double> matrix, std::size_t dimension, bool eigenvectors,
    std::optional<std::size_t> count = {}, std::string const& which = "SA");
ALPS_ED_DECL DenseEigenResult<std::complex<double>> dense_eigensystem(
    std::vector<std::complex<double>> matrix, std::size_t dimension, bool eigenvectors,
    std::optional<std::size_t> count = {}, std::string const& which = "SA");

} // namespace alps::ed
