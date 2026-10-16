// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/model/include/alps/model/qubit_encoding.hpp. Adaptation: the encoding
// snapshots a PhysicalBasis (alpsx copied SiteBasis metadata), <alps/models/...>
// include names, ALPS_MODELS_DECL exports.
//
// Qubit register layout of an ordered physical basis: every site occupies
// consecutive qubits, site 0 the least significant ones, and its local basis
// index is written as a binary or Gray code word. Unused code words of a
// non-power-of-two site are not physical states. No sector projection occurs.

#pragma once

#include <alps/models_export.h>
#include <alps/models/physical_basis.hpp>

#include <complex>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>

namespace alps {
namespace models {

class Hamiltonian;

enum class QubitEncodingScheme { binary, gray };

inline constexpr std::size_t default_local_matrix_bytes = 64 * 1024 * 1024;
inline constexpr std::size_t default_state_amplitudes = 1024 * 1024;

struct DecodedQubitState {
    std::vector<std::complex<double>> amplitudes;
    double discarded_norm_squared = 0.;
};

// Immutable snapshot of ordered local bases. Site 0 uses the least significant
// qubits; vectors of physical amplitudes use ALPS site-0-most-significant order.
// Binary/Gray code words have no amplitude phase. Finite bosonic bases use
// binary or Gray; fermionic sites must be binary occupation bases of dimension
// 2 or 4 with occupation parity, and a mixture of fermions with bosons of
// dimension above 2 is rejected (it needs a separately validated encoding).
class ALPS_MODELS_DECL QubitEncoding {
public:
    explicit QubitEncoding(PhysicalBasis const& basis,
                           QubitEncodingScheme scheme = QubitEncodingScheme::binary,
                           std::size_t max_local_matrix_bytes = default_local_matrix_bytes);
    int num_qubits() const { return num_qubits_; }
    QubitEncodingScheme scheme() const { return scheme_; }
    std::vector<int> const& site_qubits() const { return widths_; }
    std::vector<int> const& site_offsets() const { return offsets_; }
    std::vector<int> const& site_dimensions() const { return dimensions_; }
    // The snapshotted local spaces (labels, statistics, quantum numbers, parity).
    std::vector<LocalBasis> const& sites() const { return sites_; }
    std::size_t codeword(std::size_t site, int local_index) const;
    std::vector<int> encode_basis_bits(std::vector<int> const& local_indices) const;
    std::optional<std::vector<int>> decode_basis_bits(std::vector<int> const& bits) const;
    std::uint64_t state_index(std::vector<int> const& local_indices) const;
    // Throws invalid_argument unless the sites of `basis` equal the snapshot.
    void require_same_sites(PhysicalBasis const& basis) const;

    // Explicit small-state interchange, never used by operator compilation.
    // Decoding projects without renormalizing and reports discarded squared norm.
    std::vector<std::complex<double>> encode_statevector(std::vector<std::complex<double>> const& physical,
        std::size_t max_amplitudes = default_state_amplitudes) const;
    DecodedQubitState decode_statevector(std::vector<std::complex<double>> const& encoded,
        std::size_t max_amplitudes = default_state_amplitudes) const;
    // Fraction of squared norm in the physical subspace; rejects zero states.
    double physical_weight(std::vector<std::complex<double>> const& encoded,
        std::size_t max_amplitudes = default_state_amplitudes) const;

private:
    std::optional<std::size_t> physical_index(std::size_t encoded_index) const;
    std::pair<std::size_t, std::size_t> state_dimensions(std::size_t limit) const;
    std::vector<LocalBasis> sites_;
    std::vector<int> dimensions_, widths_, offsets_;
    QubitEncodingScheme scheme_;
    int num_qubits_ = 0;
};

// The encoding of physical_basis(model).
ALPS_MODELS_DECL QubitEncoding make_qubit_encoding(Hamiltonian const& model,
    QubitEncodingScheme scheme = QubitEncodingScheme::binary,
    std::size_t max_local_matrix_bytes = default_local_matrix_bytes);

} // namespace models
} // namespace alps
