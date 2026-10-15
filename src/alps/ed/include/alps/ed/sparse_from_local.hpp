// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/include/alps/diag/sparse_from_local.hpp.
// Adaptation: C++17 (no std::span), <alps/ed/...> include names, ALPS_ED_DECL exports.
//
// Sparse many-body assembly over local operators. Both the explicit-basis
// complex COO entry point and the enumerated real CSR entry point use the
// same term-application kernel. Basis enumeration stays with the caller
// when a projected Hilbert space is much smaller than its full product.
//
// Real CSR conventions (locked by tests/sparse_from_local_lock.cpp):
//   - basis states are mixed-radix codes over site_dims with SITE 0 LEAST
//     significant, enumerated ascending; a QnConstraint keeps states whose
//     per-site QN values (from diag_observables; a site type without the
//     QN contributes 0 — the kit convention) sum to `target` within 1e-9;
//   - CSR row = IN state, column = OUT state: data holds ⟨out|H|in⟩ (for
//     the real-symmetric Hamiltonians this path accepts, equal to the
//     transpose convention either way);
//   - a term's amplitude is coeff · Π_k op_k(out_k, in_k) · Π_{s∈jw} F_s(in_s)
//     with F the diagonal parity operator of the site's basis — identical
//     to the op·F / bare-F Kronecker convention of the kit's dense
//     reference assembly (src/alps/models/tests/dense_e0.hpp); a site type
//     absent from `parity_by_type` is bosonic (F = +1);
//   - out-states falling outside a constrained basis are dropped (sector
//     projection — QN-conserving Hamiltonians lose nothing);
//   - the real CSR entry point retains its existing 1e-12 imaginary-part
//     tolerance; complex COO has no such projection or tolerance.

#pragma once

#include <alps/ed_export.h>
#include <alps/models/compiled.hpp>

#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace alps {
namespace ed {

struct CooFromLocal {
    std::vector<std::int64_t> rows, cols;
    std::vector<std::complex<double>> data;
    bool is_complex = false;
};

using MatrixElementVisitor =
    std::function<void(std::int64_t row, std::int64_t col, std::complex<double> value)>;

using LocalTransitionVisitor = std::function<void(
    std::int64_t source, std::vector<std::int64_t> const& target, std::complex<double> amplitude)>;

// Apply local terms to a row-major (num_sources, site_dims.size()) batch of
// source states. Source rows may repeat. Emit uncombined nonzero transitions
// in term/source/factor-output order, including targets outside the source
// batch. This operation performs no basis lookup, projection, phase weighting
// or cutoff; those policies belong to the caller. The target vector is
// borrowed and valid only during its callback. Do not mutate input/source/
// parity data from a callback. All inputs are validated before the first
// callback; an empty visitor is invalid and callback exceptions propagate
// immediately. Returns whether any emitted amplitude or its term coefficient
// is complex.
ALPS_ED_DECL bool visit_local_transitions(
    models::CompiledOperator const& input,
    std::vector<std::int64_t> const& sources, std::size_t num_sources,
    LocalTransitionVisitor const& visitor,
    std::vector<std::vector<double>> const& parity = {});

// Stream the same nonzero, uncombined contributions as build_coo_from_local,
// retaining term order, input-column order and factor traversal order. Returns
// its is_complex flag, including complex coefficients of emitted contributions.
// The visitor must be nonempty; its exceptions propagate immediately. All input
// validation finishes before the first callback. No matrix-entry buffers are
// retained. The caller must not mutate the operator, basis or parity in a callback.
ALPS_ED_DECL bool visit_local_operator(
    models::CompiledOperator const& input,
    std::vector<std::int64_t> const& states, std::size_t num_states,
    MatrixElementVisitor const& visitor,
    std::vector<std::vector<double>> const& parity = {});

// Codes use site 0 MOST significant and must be distinct and smaller than the
// product of the dimensions, which must fit uint64_t. Codes remain in the
// caller's order. Sorted codes use binary search without an auxiliary index;
// other orders use a lookup table. Empty parity has the same meaning as above.
ALPS_ED_DECL bool visit_local_operator(
    models::CompiledOperator const& input,
    std::vector<std::uint64_t> const& codes, MatrixElementVisitor const& visitor,
    std::vector<std::vector<double>> const& parity = {});

// states is a row-major (num_states, site_dims.size()) array of distinct
// local-state indices, in any caller-selected order. No full-space
// enumeration or finite-width full-space code is required. Matrix rows
// are OUT states, columns IN states; duplicate contributions are retained
// in term order. Operators on a given site must already be composed into
// one factor. Empty parity means F=1; otherwise supply one diagonal per
// site. No imaginary component is discarded, however small.
ALPS_ED_DECL CooFromLocal build_coo_from_local(
    models::CompiledOperator const& input,
    std::vector<std::int64_t> const& states, std::size_t num_states,
    std::vector<std::vector<double>> const& parity = {});

// Keep basis states whose per-site `qn` values sum to `target`.
struct QnConstraint {
    std::string qn;
    double target = 0.0;
};

struct SparseFromLocal {
    std::vector<std::int64_t> indptr;   // dim + 1
    std::vector<std::int64_t> indices;  // nnz
    std::vector<double> data;           // nnz
    std::int64_t dim = 0;
    std::vector<std::uint64_t> states;  // mixed-radix codes, ascending
};

ALPS_ED_DECL SparseFromLocal build_sparse_from_local(
    models::CompiledModel const& model,
    std::map<int, std::vector<double>> const& parity_by_type = {},
    std::vector<QnConstraint> const& constraints = {},
    std::uint64_t max_full_dim = (std::uint64_t{1} << 24));

} // namespace ed
} // namespace alps
