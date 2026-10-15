// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/include/alps/diag/spectrum.hpp.
// Adaptation: C++17 (no std::span), <alps/ed/...> include names, ALPS_ED_DECL
// exports, namespace alps::ed; block dispatch follows the IETL-based Lanczos
// driver's policy.
#pragma once

#include <alps/ed_export.h>

#include <complex>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace alps::ed {

struct SpectrumLocalSpace {
    // One parity sign (+1 or -1) per local basis state.
    std::vector<double> parity;
    // One row per declared conserved charge, in the caller's label order.
    std::vector<std::vector<double>> quantum_numbers;
};

struct SpectrumTranslation {
    // One primitive translation, mapping source site -> translated site.
    // Generators must commute and each raised to its period must be identity.
    std::vector<std::size_t> sites;
    std::size_t period = 0;
};

struct SpectrumBlock {
    // Dimension before eigenvalue selection. Comparing this with the returned
    // eigenvalue count distinguishes a complete block from a truncated one.
    std::size_t block_dimension = 0;
    // Translation eigenphases, one per supplied generator, in [0, 2*pi).
    std::vector<double> momentum;
    std::vector<double> quantum_numbers;
    std::vector<double> eigenvalues;
    std::string method;
    // Lanczos-only effective controls/counters; absent for dense blocks.
    std::optional<double> tolerance;
    std::optional<int> max_iterations;
    std::optional<int> iterations;
    int nconv = 0;
    // Set only for finite-group blocks. All irrep rows are explicitly retained:
    // eigenvalues count physical states, not one representative per multiplet.
    std::string irrep;
    std::size_t irrep_dimension = 0;
    std::size_t group_order = 0;
};

// Resolve a Hermitian product-basis CSR matrix into joint conserved-charge
// and optional momentum blocks, including fermionic permutation signs. All
// numerical inputs are borrowed only for this call; returned blocks own their
// values. Product states are a row-major [dimension, local_spaces.size()]
// array. Charge totals use the product-basis tolerance of 1e-9. The input
// must preserve each charge block and commute with the supplied translations.
// Empty translations select charge resolution alone. Empty momentum_indices
// request every momentum; otherwise indices select rows of the Cartesian
// momentum grid (last translation direction varies fastest).
// A block is solved densely unless its dimension exceeds both 1000 and
// 8*count + 200 (the automatic eigenstate policy); SM is always dense.
// which selects SA/LA/SM/LM within every block; returned energies are ascending.
// Omitted Lanczos controls are the driver defaults, tolerance 1e-10 and
// min(10*block_dimension, 1000) steps per run. Explicit controls must be
// positive and finite, even for dense/empty blocks.
ALPS_ED_DECL std::vector<SpectrumBlock> symmetry_spectrum(
    std::vector<std::int64_t> const& indptr,
    std::vector<std::int64_t> const& indices,
    std::vector<std::complex<double>> const& values,
    std::size_t dimension,
    std::vector<std::int64_t> const& product_states,
    std::vector<SpectrumLocalSpace> const& local_spaces,
    std::vector<SpectrumTranslation> const& translations,
    int num_eigenvalues,
    std::vector<std::size_t> const& momentum_indices = {},
    std::string const& which = "SA",
    std::optional<double> tolerance = {},
    std::optional<int> max_iterations = {});

} // namespace alps::ed
