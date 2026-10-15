// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/include/alps/model/explicit_model.hpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports;
// Pauli/qubit export declarations removed.
//
// alps::models::Hamiltonian — lattice + site bases + a recipe of bond/site
// term templates. compile_hamiltonian() produces the model-owned
// CompiledModel; exact diagonalization and DMRG consume its operator.

#pragma once

#include <alps/models_export.h>
#include <alps/models/lattice.hpp>
#include <alps/models/terms.hpp>
#include <alps/models/compiled.hpp>

#include <functional>
#include <optional>

namespace alps {
namespace models {

// Shared local-term compiler for Monte Carlo frontends. Canonicalize without
// approximating local matrices, validate declared charges against the complete
// operator (including cancellations across different supports), then split
// factors into fixed-charge components. No global product basis is built.
// Charge labels are compared within 1e-9; matrix entries are never pruned by
// that tolerance. A charge must exist on at least one site; sites omitting it
// are neutral (for example localized Kondo spins carry zero particle number).
// Throws invalid_argument for an incorrect declaration.
ALPS_MODELS_DECL TermSum compile_local_terms(TermSum const& terms,
                            std::vector<SiteBasisPtr> const& sites,
                            std::vector<std::string> const& conserved_qns = {},
                            FermionConvention convention = FermionConvention::strict);

// Owned numeric form of compile_local_terms, retaining its ordered charge
// components and explicit Jordan-Wigner strings. Unlike the ED compiler below,
// this supplies the proposal-operator decomposition used by lattice samplers.
ALPS_MODELS_DECL CompiledOperator compile_local_operator(TermSum const& terms,
    std::vector<SiteBasisPtr> const& sites,
    std::vector<std::string> const& conserved_qns = {},
    FermionConvention convention = FermionConvention::strict);

// Compile a numeric Hamiltonian for model-to-solver use. The same native
// full-operator conservation check as compile_local_terms runs before any
// sector projection and like-term coefficient export. Canonical terms remain
// unsplit by QN and unpruned. Like-factor coefficients are summed exactly then
// exported once; if no finite preserving coefficient exists, retain the finite
// contributions. Descriptor term counts/order are not input-term identities.
// Same-site runs use exact intermediates followed by a double matrix export.
// Normal rounding at that export can still lose a cross-term cancellation
// remainder, including before conservation or Hermiticity validation. This
// interface does not promise exact arbitrary operator algebra.
// Charge-label tolerance is fixed at 1e-9, independently of amplitude scale.
// This does not certify Hermiticity, translation invariance or convergence;
// those remain requirements of the consuming solver.
// The raw overload fills site_dims and term/Jordan-Wigner data only: no
// lattice geometry or site-type metadata is invented. Inputs are unchanged
// and the descriptor owns its arrays. No global Hilbert space is built.
ALPS_MODELS_DECL CompiledOperator compile_hamiltonian(
    TermSum const& terms, std::vector<SiteBasisPtr> const& sites,
    std::vector<std::string> const& conserved_qns = {},
    FermionConvention convention = FermionConvention::strict);

// Lower an arbitrary observable with the operator algebra's canonicalization
// policy. Charge-changing, odd and non-Hermitian operators are allowed. This
// owns its arrays, retains exact identity semantics and applies the same local
// ordering/parity conventions as compile_hamiltonian. It imposes no eigensolver
// validation and shares the finite-export limitations documented above.
ALPS_MODELS_DECL CompiledOperator compile_operator(
    TermSum const& terms, std::vector<SiteBasisPtr> const& sites,
    FermionConvention convention = FermionConvention::strict);

// Validate the full numeric operator before a model Hermitian eigensolve.
// Does not impose charge declarations or change arbitrary operator compilation.
// Canonicalization precedes adjunction, including graded swap signs. Subject to
// the local-product export limitation above, exact accumulation recognizes
// cancellations across different site supports before like-coefficient export.
// For each independent local matrix-unit word, compare H and H^dagger using
// |a-b| <= tolerance * max(|a|,|b|). The identity word uses the largest
// nonidentity pair magnitude instead, so real energy shifts cannot hide an
// imaginary scalar. A pure scalar must be exactly real. Zero tolerance is exact.
// No global basis or Jordan-Wigner string expansion is performed. An exact
// tensor equality fast path handles long Hermitian products; the ambiguous
// sparse-word fallback throws length_error above max_expanded_terms rather
// than allocating an unbounded tensor expansion. Inputs remain unchanged.
ALPS_MODELS_DECL void validate_hermiticity(
    TermSum const& terms, std::vector<SiteBasisPtr> const& sites,
    FermionConvention convention = FermionConvention::strict,
    double tolerance = 1e-11, std::size_t max_expanded_terms = 100000);

// Bond/site term specification: a callback invoked per (source, target)
// bond or per site, returning the term to add (the Python OpSpec callable
// form; located-OpSum specs are covered by capturing lambdas).
using BondSpec = std::function<TermSum(int, int)>;
using SiteSpec = std::function<TermSum(int)>;

class ALPS_MODELS_DECL Hamiltonian {
public:
    // `sites_by_type` maps lattice site_type -> SiteBasis; the
    // single-basis overload applies one basis to every site.
    Hamiltonian(lattice::Lattice lattice, SiteBasisPtr site,
                FermionConvention conv = FermionConvention::strict,
                std::vector<std::string> conserved_qns = {});
    Hamiltonian(lattice::Lattice lattice, std::map<int, SiteBasisPtr> sites_by_type,
                FermionConvention conv = FermionConvention::strict,
                std::vector<std::string> conserved_qns = {});

    // -- term accumulation (templates expand lazily, mirroring Python) --
    Hamiltonian& add_bond(Complex coeff, BondSpec spec,
                          std::optional<int> bond_type = std::nullopt,
                          bool add_hc = false);
    Hamiltonian& add_site(Complex coeff, SiteSpec spec,
                          std::optional<int> site_type = std::nullopt,
                          bool add_hc = false);
    Hamiltonian& add(TermSum term, bool add_hc = false);

    // -- introspection --
    lattice::Lattice const& lat() const { return lattice_; }
    int num_sites() const { return lattice_.num_sites(); }
    SiteBasisPtr site_for(int lattice_index) const;
    std::vector<int> site_dims() const;
    std::vector<std::string> const& conserved_qns() const { return conserved_qns_; }
    FermionConvention fermion_convention() const { return convention_; }

    // Expand every template into one TermSum (not canonicalised) —
    // template-major, bond/site-minor, in insertion order.
    TermSum expand() const;

    // The Monte-Carlo descriptor: expand -> compile_local_terms (including
    // full-operator conservation validation) -> per term, emit the
    // canonical factors, effective complex coefficient, and the explicit
    // Jordan-Wigner F-string sites; plus site_dims/site_types, the lattice
    // edges, and per-site-type diag observables (each declared QN's
    // per-state values — what alps::qmc::DiagonalObservables::build_diagonal_operator consumes).
    CompiledModel to_local_decomposition() const;

    // Unfiltered local operator, with this model's geometry and fermion
    // convention. These preserve every charge component without requiring
    // the Hamiltonian's symmetry declaration. Use them for arbitrary (including
    // charge-changing or non-Hermitian) observables. The no-argument overload
    // lowers the full Hamiltonian. No product basis is constructed.
    CompiledModel operator_decomposition() const;
    CompiledModel operator_decomposition(TermSum const& op) const;

private:
    CompiledModel lower_operator(TermSum const& canonical, bool jordan_wigner) const;
    struct BondTemplate { Complex coeff; BondSpec spec; std::optional<int> bond_type; bool add_hc; };
    struct SiteTemplate { Complex coeff; SiteSpec spec; std::optional<int> site_type; bool add_hc; };
    struct GlobalTerm  { TermSum op; bool add_hc; };

    lattice::Lattice lattice_;
    SiteBasisPtr homogeneous_site_;                 // nullptr when per-type map
    std::map<int, SiteBasisPtr> site_map_;
    FermionConvention convention_;
    std::vector<std::string> conserved_qns_;
    std::vector<BondTemplate> bond_templates_;
    std::vector<SiteTemplate> site_templates_;
    std::vector<GlobalTerm> global_terms_;
};

// Expand and compile the model under its declared charges and fermion
// convention, including its lattice graph metadata and per-site-type diagonal
// metadata. Arbitrary charge-changing observables continue to use
// Hamiltonian::operator_decomposition instead.
ALPS_MODELS_DECL CompiledModel compile_hamiltonian(Hamiltonian const& model);

ALPS_MODELS_DECL void validate_hermiticity(Hamiltonian const& model,
    double tolerance = 1e-11, std::size_t max_expanded_terms = 100000);

// The model's ordered physical local spaces, declared charges and convention.
ALPS_MODELS_DECL PhysicalBasis physical_basis(Hamiltonian const& model);

// Sites needing a Jordan-Wigner F for canonical-ordered factors: every site
// k in [0, n_sites) where the count of odd-parity factor sites strictly
// greater than k is odd (compile._jw_chain_sites).
ALPS_MODELS_DECL std::vector<int> jw_chain_sites(std::vector<Factor> const& factors, int n_sites);

} // namespace models
} // namespace alps
