// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, the PauliOperator,
// PauliExportOptions and export declarations of
// src/model/include/alps/model/explicit_model.hpp. Adaptation: the export
// takes a CompiledOperator on its PhysicalBasis (alpsx took raw terms with
// their site bases and canonicalised them itself); charges, convention and
// local spaces come from the basis; a require_conservation switch replaces
// the empty-charge-list idiom for observables; <alps/models/...> include
// names, ALPS_MODELS_DECL exports.
//
// Sparse Pauli expansion of a compiled operator in the qubit register of a
// QubitEncoding. Strings put qubit 0 at the RIGHT end. Finite bosonic bases
// use binary/Gray code words; strict dim-2/4 fermions retain binary occupation
// parity and their Jordan-Wigner strings become Z strings. Unused code words
// of a non-power-of-two site carry the trace-completed local identity. No
// sector projection or chemistry-integral conversion.

#pragma once

#include <alps/models_export.h>
#include <alps/models/compiled.hpp>
#include <alps/models/physical_basis.hpp>
#include <alps/models/qubit_encoding.hpp>

#include <complex>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace alps {
namespace models {

class Hamiltonian;

struct ALPS_MODELS_DECL PauliOperator {
    int num_qubits = 0;
    std::vector<int> site_qubits;
    // Site-major CAR order; spinful sites contribute (dn, up), i.e. the
    // higher local bit first. Empty for bosonic sites.
    std::vector<int> fermion_mode_qubits;
    std::map<std::string, std::complex<double>> terms;
    // Sum of absolute discarded Pauli coefficients, including imaginary
    // roundoff removed when Hermiticity was requested (operator-norm bound).
    double discarded_l1 = 0.;
    std::optional<QubitEncoding> encoding;
    std::size_t local_decompositions = 0;
    std::size_t local_cache_hits = 0;
    // Register index of a physical product state (one local index per site).
    std::uint64_t state_index(std::vector<int> const& local_indices) const;
    double coefficient_l1(bool include_identity = true) const;
};

struct PauliExportOptions {
    double tolerance = 1e-12;
    // Bounds both each term's Cartesian expansion and the accumulated map.
    std::size_t max_terms = 100000;
    bool require_hermitian = true;
    bool require_even = true;
    // Check that the operator commutes with every conserved charge of the
    // basis. Switch the three requirements off for an arbitrary observable.
    bool require_conservation = true;
    // Bounds each dense local transform; the per-call cache has the same
    // byte budget (estimated including map nodes), and at most 64 entries.
    std::size_t max_local_matrix_bytes = default_local_matrix_bytes;
    // Separately control Hermiticity/parity/charge roundoff checks when zero
    // coefficient pruning is wanted. Omitted retains the pruning tolerance.
    std::optional<double> validation_tolerance;
};

// The operator's local matrices and Jordan-Wigner strings in the register of
// `encoding` (or of the binary encoding of `basis`). Every term may act on any
// number of sites; a site repeated within a term multiplies its factors in
// order. Fermionic sites need the strict convention. Throws invalid_argument
// for a basis/encoding/operator mismatch or a failed requirement and
// length_error above a resource limit.
ALPS_MODELS_DECL PauliOperator export_pauli(CompiledOperator const& op, PhysicalBasis const& basis,
                                            PauliExportOptions const& options = {});
ALPS_MODELS_DECL PauliOperator export_pauli(CompiledOperator const& op, PhysicalBasis const& basis,
                                            QubitEncoding const& encoding,
                                            PauliExportOptions const& options = {});
// compile_hamiltonian(model).op on physical_basis(model).
ALPS_MODELS_DECL PauliOperator export_pauli(Hamiltonian const& model, PauliExportOptions const& options = {});
ALPS_MODELS_DECL PauliOperator export_pauli(Hamiltonian const& model, QubitEncoding const& encoding,
                                            PauliExportOptions const& options = {});

// Projectors are constructed directly in the register, not by extending the
// physical identity (whose identity-preserving extension would be I).
ALPS_MODELS_DECL PauliOperator local_valid_projector(QubitEncoding const& encoding, std::size_t site,
                                                     PauliExportOptions const& options = {});
// Explicit lambda * sum_i (I-P_i); does not modify an exported Hamiltonian.
ALPS_MODELS_DECL PauliOperator invalid_state_penalty(QubitEncoding const& encoding, double strength,
                                                     PauliExportOptions const& options = {});

} // namespace models
} // namespace alps
