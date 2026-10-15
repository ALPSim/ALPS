// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/include/alps/model/site.hpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.
//
// Local Hilbert spaces shared by the scientific solver interfaces.
// Independent matrix references are in tests/models_kit_parity.cpp.
//
// A SiteBasis declares the local dimension, exchange statistics, the named
// LocalOp registry (dense matrices + fermionic parity tags + per-QN deltas),
// per-state quantum numbers, and the diagonal (-1)^F parity operator.

#pragma once

#include <alps/models_export.h>
#include <alps/models/physical_basis.hpp>

#include <complex>
#include <cstddef>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace alps {
namespace models {

using Complex = std::complex<double>;
// Row-major dim x dim dense matrix.
using Matrix = std::vector<Complex>;

enum class Parity { even, odd };

// Infer per-QN deltas of a matrix: for every declared QN q, if all non-zero
// entries M[r,c] share one q[r]-q[c] (within tol), the result maps q -> that
// delta; disagreeing QNs are omitted. The zero matrix respects every QN with
// delta 0. tol compares charge labels only: every nonzero matrix entry
// participates, however small, since an external coefficient may amplify it.
ALPS_MODELS_DECL std::map<std::string, double> compute_qn_deltas(
    Matrix const& matrix, int dim,
    std::map<std::string, std::vector<double>> const& quantum_numbers,
    double tol = 1e-9);

// A named local operator.
struct ALPS_MODELS_DECL LocalOp {
    std::string name;
    int dim = 0;
    Matrix matrix;                       // row-major dim x dim
    Parity parity = Parity::even;
    std::map<std::string, double> qn_deltas;

    Complex operator()(int r, int c) const {
        return matrix[static_cast<std::size_t>(r) * dim + c];
    }
    bool is_identity(double tol = 1e-12) const;

    // Hermitian conjugate: M -> M^dagger, deltas negate, name gains/loses
    // the "_dag" suffix.
    LocalOp dag() const;

    // Decompose into QN-delta-fixed components (one per distinct delta
    // tuple over the declared QNs, in sorted delta-tuple order; components
    // sum to the original matrix). Components use exact differences of the
    // supplied double charge labels, never rounded/tolerantly inferred labels.
    // qn_deltas and "name[Sz=+1]" tags are floating-point presentation metadata.
    std::vector<LocalOp> split_by_qn_delta(
        std::map<std::string, std::vector<double>> const& quantum_numbers) const;
};

// A local Hilbert space. NOTE: quantum_numbers is a sorted map, whereas the
// Python spec preserves declaration order. Every built-in site declares its
// QNs in alphabetical order (Sz; N; N,Sz), so ordering-observable behavior
// (split bucket tags, delta tuples) agrees; custom sites with non-sorted QN
// declaration order are the one place the two implementations may order
// split components differently (the component *set* is identical).
struct ALPS_MODELS_DECL SiteBasis {
    int dim = 0;
    Statistics statistics = Statistics::bosonic;
    std::map<std::string, LocalOp> operators;
    std::map<std::string, std::vector<double>> quantum_numbers;
    std::vector<double> parity_operator;   // diagonal (-1)^F, length dim
    // Optional ordered physical-state labels for interchange. When absent,
    // qubit encodings use ordinal labels; quantum numbers still define order.
    std::vector<std::string> basis_labels;

    LocalOp const& op(std::string const& name) const;   // throws on miss
    bool has_op(std::string const& name) const { return operators.count(name) != 0; }
};

// Shared-ownership handle: Term factors reference their SiteBasis.
using SiteBasisPtr = std::shared_ptr<SiteBasis const>;

// Finalise a hand-built SiteBasis exactly like Site.__init__: add the
// identity operator "I", default the parity operator to +1, and back-fill
// missing qn_deltas on every operator from the declared quantum numbers.
ALPS_MODELS_DECL SiteBasisPtr finalize_site(SiteBasis site);

// Resolve ordered local spaces into a complete, owned physical basis. Each
// conserved charge must exist on at least one site; a site that omits such a
// charge is neutral and receives an explicit zero row. Other quantum numbers,
// labels, statistics and parity are preserved unchanged. Missing optional
// SiteBasis labels become ordinal strings ("0", "1", ...). No operator
// matrices, product basis or solver runtime are required. Inputs are not
// modified. Complete PhysicalBasis values and archives remain subject to the
// stricter validate_basis contract; they are never repaired here.
ALPS_MODELS_DECL PhysicalBasis physical_basis(
    std::vector<SiteBasisPtr> const& sites,
    std::vector<std::string> const& conserved_qns = {},
    FermionConvention convention = FermionConvention::strict);

// Check model-site metadata against an existing complete physical basis.
// The reference is validated once for the entire batch. Only locally missing
// declared conserved rows may be supplied as zeros; all other metadata must
// match exactly. The reference already establishes the charge context, so a
// subset need not contain a site on which every charge was originally declared.
// This does not change complete-basis equality or archive validation. Neither
// the reference nor the supplied model sites are modified.
ALPS_MODELS_DECL void require_compatible_sites(
    PhysicalBasis const& reference,
    std::vector<std::pair<std::size_t, SiteBasisPtr>> const& sites);

// -- built-in sites (mirroring sites.py) --------------------------------------
ALPS_MODELS_DECL SiteBasisPtr spin(double S = 0.5);           // Sz,Sp,Sm,Sx,Sy,S2 ; QN Sz
ALPS_MODELS_DECL SiteBasisPtr boson(int n_max = 3);           // a,a_dag,n (+b aliases) ; QN N
ALPS_MODELS_DECL SiteBasisPtr hardcore_boson();               // a,a_dag,n (+b aliases) ; QN N
ALPS_MODELS_DECL SiteBasisPtr spinless_fermion();             // c,c_dag,n,F ; QN N
ALPS_MODELS_DECL SiteBasisPtr spinful_fermion();              // c_*,n_*,S* dim 4 ; QN N,Sz

} // namespace models
} // namespace alps
