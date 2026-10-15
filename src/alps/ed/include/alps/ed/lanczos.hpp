// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Interface after skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/include/alps/diag/lanczos.hpp,
// whose implementation used Spectra. This implementation is new: a deflated
// Lanczos driver over the in-tree IETL (src/alps/ietl), see src/lanczos.cpp.
//
// Contract:
//   - real symmetric or complex Hermitian CSR; duplicate and unordered
//     entries are accumulated; offsets, indices and values are validated;
//   - 0 < count < dimension; `which` is SA (lowest), LA (highest) or LM
//     (largest magnitude). SM is not available: Lanczos without shift-invert
//     does not resolve interior eigenvalues; use dense_eigensystem;
//   - tolerance <= 0 selects 1e-10, max_iterations <= 0 selects
//     min(10*dimension, 1000) Lanczos steps per run (never fewer than 2);
//   - degenerate eigenvalues are resolved by locking converged eigenvectors
//     and restarting on the deflated operator, so `count` eigenpairs come
//     back with multiplicity, orthonormal to roundoff;
//   - eigenvalues are ascending; eigenvectors are column-major
//     (eigenvectors[column * dimension + row]) and only present when
//     requested; residual_norms holds ||A v - lambda v||_2 per column;
//   - non-convergence or a Lanczos breakdown throws std::runtime_error. The
//     method is meant for dimensions well above `count`; small spaces belong
//     to dense_eigensystem.

#pragma once

#include <alps/ed_export.h>

#include <complex>
#include <cstdint>
#include <string>
#include <vector>

namespace alps::ed {

template<class Scalar>
struct LanczosResult {
    std::vector<double> eigenvalues;        // ascending, `count` values
    std::vector<Scalar> eigenvectors;       // column-major dimension x count, when requested
    std::vector<double> residual_norms;     // per column, when vectors were requested
    int nconv = 0;
    double tolerance = 0.;
    int max_iterations = 0;                 // per-run Lanczos step limit in effect
    int iterations = 0;                     // Lanczos steps summed over all runs
    int runs = 0;                           // Lanczos runs (initial plus deflated restarts)
};

ALPS_ED_DECL LanczosResult<double> lanczos_csr(
    std::vector<std::int64_t> const& offsets, std::vector<std::int64_t> const& indices,
    std::vector<double> const& values, std::int64_t dimension, int count,
    std::string const& which = "SA", double tolerance = 0., int max_iterations = 0,
    bool want_vectors = true);
ALPS_ED_DECL LanczosResult<std::complex<double>> lanczos_csr(
    std::vector<std::int64_t> const& offsets, std::vector<std::int64_t> const& indices,
    std::vector<std::complex<double>> const& values, std::int64_t dimension, int count,
    std::string const& which = "SA", double tolerance = 0., int max_iterations = 0,
    bool want_vectors = true);

} // namespace alps::ed
