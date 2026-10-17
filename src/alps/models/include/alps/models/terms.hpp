// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/include/alps/model/terms.hpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.
//
// Native operator algebra. A Term is one coefficient times an ordered product of located local
// operators with an accumulated fermionic swap sign; a TermSum is a formal
// linear combination with distributive arithmetic, Hermitian conjugation,
// canonicalisation, and QN-delta analysis.
//
// Coefficients are bound std::complex<double> values. Python resolves its
// symbolic parameters before calling the native compiler.

#pragma once

#include <alps/models_export.h>
#include <alps/models/site.hpp>

#include <vector>

namespace alps {
namespace models {

// One product factor: a local operator located at a lattice site index.
// Carries its SiteBasis so canonicalisation can merge same-site factors
// (re-inferring QN deltas) and the QN split can consult declared QNs.
struct Factor {
    int site = 0;
    LocalOp op;
    SiteBasisPtr basis;
};

// coeff * (op@i1 * op@i2 * ...) with sign in {+1,-1} accumulated from
// cross-site odd-odd reordering. Effective coefficient = coeff * sign.
struct ALPS_MODELS_DECL Term {
    Complex coeff{1.0, 0.0};
    std::vector<Factor> factors;
    int sign = 1;

    Complex effective_coeff() const { return sign == 1 ? coeff : -coeff; }

    // Per-QN total delta: intersection of QN names every factor respects,
    // values summed. Empty product -> empty map.
    std::map<std::string, double> qn_delta() const;
};

// Reorder factors to non-decreasing site index (stable bubble sort; each
// cross-site swap of two odd-parity operators flips the sign when
// apply_swap_sign — the fermion anticommutation rule), then merge same-site
// runs into single matrix-product operators (parity XORs, QN deltas
// re-inferred from the parent site), dropping identity factors. Complete
// same-site products retain exact finite-double intermediates until their
// final double export. Powers of two may move between matrices and coefficient
// to avoid range loss; a product without such a preserving scale throws.
// Normal rounding at final export remains and can affect later cancellation.
ALPS_MODELS_DECL Term canonicalise_term(Term const& t, bool apply_swap_sign = true);
// Explicit identity tolerance for compilers that account for all pruning
// themselves. The two-argument overload retains its historical 1e-12 default.
ALPS_MODELS_DECL Term canonicalise_term(Term const& t, bool apply_swap_sign, double identity_tolerance);

class ALPS_MODELS_DECL TermSum {
public:
    TermSum() = default;
    explicit TermSum(std::vector<Term> terms) : terms_(std::move(terms)) {}

    static TermSum zero() { return TermSum(); }
    // A located single-operator term.
    static TermSum located(SiteBasisPtr basis, std::string const& op_name, int site);

    std::vector<Term> const& terms() const { return terms_; }
    bool empty() const { return terms_.empty(); }

    // -- algebra (mirrors OpSum arithmetic) --
    TermSum operator+(TermSum const& other) const;
    // In-place sums append the other operand's terms, in order, without copying
    // this sum's terms; accumulating N terms is linear. `s += s` is safe.
    TermSum& operator+=(TermSum const& other);
    TermSum& operator+=(TermSum&& other);
    TermSum& operator-=(TermSum const& other);
    TermSum operator-(TermSum const& other) const;
    TermSum operator-() const;
    TermSum operator*(TermSum const& other) const;   // distributive concat
    TermSum operator*(Complex scalar) const;
    friend TermSum operator*(Complex scalar, TermSum const& s) { return s * scalar; }

    // Hermitian conjugate: per term, reverse factor order + dagger each
    // local op + conjugate the effective coefficient (sign folded to +1).
    TermSum dag() const;
    TermSum add_hc() const { return *this + dag(); }

    // Canonicalise every term and combine like terms (keyed on the factor
    // structure: site, numeric matrix entries (signed zeros equal), parity;
    // first-appearance group order kept; exactly-zero coefficients dropped).
    // Like coefficients accumulate exactly before one double export. If that
    // export would overflow or lose a subnormal component, retain the original
    // finite contributions instead. Normal final rounding remains possible.
    TermSum canonicalise(bool apply_swap_sign = true) const;
    // Faithful lowering uses zero identity tolerance before applying external
    // coefficients; the one-argument overload retains its existing policy.
    TermSum canonicalise(bool apply_swap_sign, double identity_tolerance) const;

    // Decompose every factor into QN-delta-fixed components (Cartesian
    // product across factors, component order = sorted delta tuples).
    TermSum split_factors_by_qn_delta() const;

private:
    std::vector<Term> terms_;
};

} // namespace models
} // namespace alps
