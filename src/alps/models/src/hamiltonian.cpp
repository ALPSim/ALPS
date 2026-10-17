// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/src/operators/hamiltonian.cpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.
//
// Implementation of alps::models::Hamiltonian. Mirrors hamiltonian.py
// (template expansion) and compile.py (to_local_decomposition); parity
// pinned by tests/models_kit_parity.cpp.

#include <alps/models/hamiltonian.hpp>
#include "exact_arithmetic.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <set>
#include <stdexcept>

namespace alps {
namespace models {

Hamiltonian::Hamiltonian(lattice::Lattice lattice, SiteBasisPtr site,
                         FermionConvention conv,
                         std::vector<std::string> conserved_qns)
    : lattice_(std::move(lattice)), homogeneous_site_(std::move(site)),
      convention_(conv), conserved_qns_(std::move(conserved_qns)) {
    if (!homogeneous_site_)
        throw std::invalid_argument("Hamiltonian: null SiteBasis");
}

Hamiltonian::Hamiltonian(lattice::Lattice lattice,
                         std::map<int, SiteBasisPtr> sites_by_type,
                         FermionConvention conv,
                         std::vector<std::string> conserved_qns)
    : lattice_(std::move(lattice)), site_map_(std::move(sites_by_type)),
      convention_(conv), conserved_qns_(std::move(conserved_qns)) {
    for (auto const& s : lattice_.sites())
        if (!site_map_.count(s.site_type))
            throw std::invalid_argument(
                "Hamiltonian: sites map missing entry for site_type " +
                std::to_string(s.site_type));
}

Hamiltonian& Hamiltonian::add_bond(Complex coeff, BondSpec spec,
                                   std::optional<int> bond_type, bool add_hc) {
    bond_templates_.push_back({coeff, std::move(spec), bond_type, add_hc});
    return *this;
}

Hamiltonian& Hamiltonian::add_site(Complex coeff, SiteSpec spec,
                                   std::optional<int> site_type, bool add_hc) {
    site_templates_.push_back({coeff, std::move(spec), site_type, add_hc});
    return *this;
}

Hamiltonian& Hamiltonian::add(TermSum term, bool add_hc) {
    global_terms_.push_back({std::move(term), add_hc});
    return *this;
}

SiteBasisPtr Hamiltonian::site_for(int lattice_index) const {
    if (homogeneous_site_) return homogeneous_site_;
    int const stype =
        lattice_.sites()[static_cast<std::size_t>(lattice_index)].site_type;
    return site_map_.at(stype);
}

std::vector<int> Hamiltonian::site_dims() const {
    std::vector<int> out;
    out.reserve(static_cast<std::size_t>(num_sites()));
    for (int i = 0; i < num_sites(); ++i) out.push_back(site_for(i)->dim);
    return out;
}

// Accumulated in place: rebuilding the sum per added term (total = total + t)
// copied every earlier term each time, quadratic in the term count.
TermSum Hamiltonian::expand() const {
    TermSum total;
    for (BondTemplate const& tpl : bond_templates_) {
        for (auto const& bond : lattice_.bonds()) {
            if (tpl.bond_type && bond.bond_type != *tpl.bond_type) continue;
            TermSum contribution = tpl.spec(bond.source, bond.target) * tpl.coeff;
            if (tpl.add_hc) contribution += contribution.dag();
            total += std::move(contribution);
        }
    }
    for (SiteTemplate const& tpl : site_templates_) {
        for (auto const& site : lattice_.sites()) {
            if (tpl.site_type && site.site_type != *tpl.site_type) continue;
            TermSum contribution = tpl.spec(site.index) * tpl.coeff;
            if (tpl.add_hc) contribution += contribution.dag();
            total += std::move(contribution);
        }
    }
    for (GlobalTerm const& g : global_terms_) {
        total += g.op;
        if (g.add_hc) total += g.op.dag();
    }
    return total;
}

std::vector<int> jw_chain_sites(std::vector<Factor> const& factors, int n_sites) {
    std::vector<int> odd_sites;
    for (Factor const& f : factors)
        if (f.op.parity == Parity::odd) odd_sites.push_back(f.site);
    std::vector<int> out;
    if (odd_sites.empty()) return out;
    for (int k = 0; k < n_sites; ++k) {
        int n_after = 0;
        for (int s : odd_sites)
            if (s > k) ++n_after;
        if (n_after % 2 == 1) out.push_back(k);
    }
    return out;   // ascending by construction
}

namespace {
constexpr double charge_tolerance = 1e-9;

bool finite(Complex value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
}

void validate_local_inputs(TermSum const& terms, std::vector<SiteBasisPtr> const& sites,
                           std::vector<std::string> const& conserved) {
    for (auto const& site : sites) {
        if (!site || site->dim <= 0)
            throw std::invalid_argument("local compilation requires positive site dimensions");
        if (site->parity_operator.size() != std::size_t(site->dim))
            throw std::invalid_argument("local compilation parity dimension mismatch");
        for (auto value : site->parity_operator)
            if (value != 1. && value != -1.)
                throw std::invalid_argument("local compilation parity values must be +1 or -1");
        for (auto const& [name, values] : site->quantum_numbers) {
            if (values.size() != std::size_t(site->dim))
                throw std::invalid_argument("local compilation quantum-number dimension mismatch");
            for (double value : values)
                if (!std::isfinite(value))
                    throw std::invalid_argument("local compilation requires finite quantum numbers");
        }
    }
    for (auto const& name : conserved)
        if (std::none_of(sites.begin(), sites.end(), [&](auto const& site) {
                return (site->quantum_numbers.count(name) != 0);
            }))
            throw std::invalid_argument("no site defines declared conserved quantum number " + name);
    for (auto const& term : terms.terms()) {
        if (!finite(term.coeff) || (term.sign != 1 && term.sign != -1))
            throw std::invalid_argument("local compilation requires finite coefficients and valid signs");
        for (auto const& factor : term.factors) {
            if (factor.site < 0 || std::size_t(factor.site) >= sites.size())
                throw std::invalid_argument("local operator site is out of range");
            auto const& site = *sites[std::size_t(factor.site)];
            if (factor.op.dim != site.dim
                || factor.op.matrix.size() != std::size_t(site.dim) * site.dim)
                throw std::invalid_argument("local operator matrix dimension mismatch");
            if (!factor.basis || factor.basis->dim != site.dim
                || factor.basis->statistics != site.statistics
                || factor.basis->quantum_numbers != site.quantum_numbers
                || factor.basis->parity_operator != site.parity_operator
                || factor.basis->basis_labels != site.basis_labels)
                throw std::invalid_argument("local operator basis mismatch");
            for (int row = 0; row < site.dim; ++row) for (int column = 0; column < site.dim; ++column) {
                auto value = factor.op.matrix[std::size_t(row) * site.dim + column];
                if (!finite(value))
                    throw std::invalid_argument("local compilation requires finite matrix entries");
                double const parity = factor.op.parity == Parity::odd ? -1. : 1.;
                if (value != Complex{}
                    && site.parity_operator[row] * site.parity_operator[column] != parity)
                    throw std::invalid_argument("local matrix disagrees with its fermionic parity tag");
            }
        }
    }
}

// A tensor word in the independent local basis {I, E_00,...,E_(d-2,d-2),
// E_rc (r != c)}. Identities are omitted from the word. This representation
// makes E_00 + E_11 = I exact for a two-state site, so terms on different
// supports can cancel without extending them to a global dense matrix.
using Word = std::vector<std::array<int, 3>>; // (site, row, column)
using detail::Exact;
using detail::ExactComplex;

bool changes_charge(Exact const& delta) {
    Exact const tolerance(charge_tolerance);
    return delta > tolerance || delta < -tolerance;
}

bool changes_charge(Word const& word, std::vector<SiteBasisPtr> const& sites,
                    std::string const& name) {
    Exact delta = 0;
    for (auto const& [site, row, column] : word) {
        auto const& numbers = sites[std::size_t(site)]->quantum_numbers;
        auto found = numbers.find(name);
        if (found == numbers.end()) continue; // Neutral local space, e.g. Kondo spins and N.
        auto const& charge = found->second;
        delta += Exact(charge[std::size_t(row)]) - Exact(charge[std::size_t(column)]);
    }
    return changes_charge(delta);
}

bool individually_conserving(Term const& term, std::vector<std::string> const& conserved) {
    // Common physical charges are modest integers/half-integers. Certify
    // those with integer arithmetic; every other label takes the exact sparse
    // word path. Tolerant per-factor inference cannot certify a total delta:
    // several small local spreads may add to more than charge_tolerance.
    constexpr std::int64_t bound = std::int64_t{1} << 52;
    auto twice_charge = [](double value) -> std::optional<std::int64_t> {
        double const twice = 2. * value;
        if (std::abs(twice) > double(std::int64_t{1} << 51) || std::trunc(twice) != twice)
            return std::nullopt;
        return static_cast<std::int64_t>(twice);
    };
    for (auto const& name : conserved) {
        std::int64_t total = 0;
        for (auto const& factor : term.factors) {
            auto found = factor.basis->quantum_numbers.find(name);
            if (found == factor.basis->quantum_numbers.end()) continue;
            auto const& charges = found->second;
            std::optional<std::int64_t> delta;
            for (int row = 0; row < factor.op.dim; ++row)
                for (int column = 0; column < factor.op.dim; ++column) {
                    if (factor.op(row, column) == Complex{}) continue;
                    auto const left = twice_charge(charges[std::size_t(row)]);
                    auto const right = twice_charge(charges[std::size_t(column)]);
                    if (!left || !right) return false;
                    auto const current = *left - *right;
                    if (delta && current != *delta) return false;
                    delta = current;
                }
            if (!delta) return true; // A zero local factor makes the term zero.
            total += *delta;
            if (total < -bound || total > bound) return false;
        }
        if (total != 0) return false;
    }
    return true;
}

void validate_full_conservation(TermSum const& terms,
                                std::vector<SiteBasisPtr> const& sites,
                                std::vector<std::string> const& conserved) {
    if (conserved.empty()) return;
    std::map<Word, ExactComplex> residuals;
    for (auto const& term : terms.terms()) {
        if (!finite(term.coeff))
            throw std::invalid_argument("non-finite combined local coefficient");
        if (individually_conserving(term, conserved)) continue;

        std::map<int, Matrix> local;
        for (auto const& factor : term.factors) local[factor.site] = factor.op.matrix;
        // Checked homogeneous parity makes each matrix-unit word determine
        // its odd sites uniquely. All contributions to that word therefore
        // share the same invertible JW string, which commutes with diagonal
        // charges and factors out of the zero test. Expanding that string
        // would spuriously make a long-range two-site check exponential.
        std::vector<std::pair<Word, ExactComplex>> expansion{{{}, ExactComplex(term.effective_coeff())}};
        for (auto const& [index, matrix] : local) {
            int const dim = sites[std::size_t(index)]->dim;
            ExactComplex const identity(matrix[std::size_t(dim) * dim - 1]);
            std::vector<std::pair<Word, ExactComplex>> next;
            for (auto const& [word, coefficient] : expansion) {
                if (!identity.zero()) next.emplace_back(word, coefficient * identity);
                for (int row = 0; row < dim; ++row) for (int column = 0; column < dim; ++column) {
                    ExactComplex value(matrix[std::size_t(row) * dim + column]);
                    if (row == column) value = value - identity;
                    if (value.zero()) continue;
                    Word extended = word;
                    extended.push_back({index, row, column});
                    next.emplace_back(std::move(extended), coefficient * value);
                }
            }
            expansion = std::move(next);
        }
        for (auto const& [word, value] : expansion) {
            bool const changes = std::any_of(conserved.begin(), conserved.end(), [&](auto const& name) {
                return changes_charge(word, sites, name);
            });
            if (!changes) continue;
            residuals[word] += value;
        }
    }
    for (auto const& [word, sum] : residuals) {
        if (sum.zero()) continue;
        for (auto const& name : conserved)
            if (changes_charge(word, sites, name))
                throw std::invalid_argument("Hamiltonian does not conserve declared quantum number " + name);
    }
}

TermSum validated_hamiltonian_terms(TermSum const& terms, std::vector<SiteBasisPtr> const& sites,
                                    std::vector<std::string> const& conserved_qns,
                                    FermionConvention convention) {
    validate_local_inputs(terms, sites, conserved_qns);
    bool const strict = convention == FermionConvention::strict;
    // Validate before the final double export of like-term coefficients.
    // Same-site runs retain exact intermediates, with one local-matrix export.
    std::vector<Term> canonical;
    canonical.reserve(terms.terms().size());
    for (auto const& term : terms.terms()) canonical.push_back(canonicalise_term(term, strict, 0.));
    TermSum unmerged(std::move(canonical));
    validate_local_inputs(unmerged, sites, conserved_qns);
    validate_full_conservation(unmerged, sites, conserved_qns);
    return unmerged;
}

// std::erase_if for maps is C++20; the SDK is C++17.
template <class Map>
void erase_zero_entries(Map& map) {
    for (auto it = map.begin(); it != map.end();)
        it = it->second.zero() ? map.erase(it) : std::next(it);
}

using TensorKey = std::vector<std::pair<int, std::vector<double>>>;
std::map<TensorKey, ExactComplex> exact_tensors(TermSum const& terms) {
    std::map<TensorKey, ExactComplex> result;
    for (auto const& term : terms.terms()) {
        TensorKey key;
        for (auto const& factor : term.factors) {
            std::vector<double> values;
            values.reserve(2 * factor.op.matrix.size());
            for (auto z : factor.op.matrix) {
                values.push_back(z.real());
                values.push_back(z.imag());
            }
            // Numeric keys equate signed zeros introduced by conjugation.
            key.emplace_back(factor.site, std::move(values));
        }
        result[std::move(key)] += ExactComplex(term.effective_coeff());
    }
    erase_zero_entries(result);
    return result;
}

std::map<Word, ExactComplex> hermitian_words(TermSum const& terms, std::size_t limit) {
    std::map<Word, ExactComplex> result;
    auto check_size = [limit](std::size_t size) {
        if (size > limit)
            throw std::length_error("Hermiticity validation exceeds max_expanded_terms");
    };
    for (auto const& term : terms.terms()) {
        if (term.effective_coeff() == Complex{}) continue;
        std::vector<std::pair<Word, ExactComplex>> expansion{{{}, ExactComplex(term.effective_coeff())}};
        for (auto const& factor : term.factors) {
            int const dim = factor.op.dim;
            ExactComplex const identity(factor.op.matrix.back());
            std::vector<std::pair<Word, ExactComplex>> next;
            for (auto const& [word, coefficient] : expansion) {
                if (!identity.zero()) {
                    check_size(next.size() + 1);
                    next.emplace_back(word, coefficient * identity);
                }
                for (int row = 0; row < dim; ++row) for (int column = 0; column < dim; ++column) {
                    ExactComplex value(factor.op(row, column));
                    if (row == column) value = value - identity;
                    if (value.zero()) continue;
                    check_size(next.size() + 1);
                    Word extended = word;
                    extended.push_back({factor.site, row, column});
                    next.emplace_back(std::move(extended), coefficient * value);
                }
            }
            expansion = std::move(next);
        }
        for (auto const& [word, value] : expansion) {
            auto [found, inserted] = result.try_emplace(word);
            if (inserted) check_size(result.size());
            found->second += value;
        }
    }
    erase_zero_entries(result);
    return result;
}

CompiledOperator emit_operator(TermSum const& full, int num_sites, bool jordan_wigner) {
    CompiledOperator out;
    out.terms.reserve(full.terms().size());
    for (Term const& raw : full.terms()) {
        // Preserve canonicalization, coefficient association and factor order.
        Term const t = canonicalise_term(raw, /*apply_swap_sign=*/true, 0.);
        LocalTerm term;
        term.coefficient = (t.sign == 1 ? Complex(1.0) : Complex(-1.0)) * t.coeff;
        term.factors.reserve(t.factors.size());
        for (Factor const& factor : t.factors)
            term.factors.push_back({factor.site, factor.op.matrix});
        if (jordan_wigner) term.jordan_wigner = jw_chain_sites(t.factors, num_sites);
        out.terms.push_back(std::move(term));
    }
    return out;
}

void add_model_metadata(CompiledModel& out, Hamiltonian const& model) {
    out.op.site_dims = model.site_dims();
    out.site_types = model.lat().site_types();
    out.edges = model.lat().edges();

    // Generic model metadata contains QN labels from the first site of each
    // type. QMC preparation additionally resolves named operators and physical
    // bases in its domain-level prepare_model contract.
    std::set<int> seen;
    for (int i = 0; i < model.num_sites(); ++i) {
        int const stype = out.site_types[static_cast<std::size_t>(i)];
        if (!seen.insert(stype).second) continue;
        SiteBasisPtr const site = model.site_for(i);
        if (site->quantum_numbers.empty()) continue;
        auto& dst = out.diag_observables[stype];
        for (auto const& [name, arr] : site->quantum_numbers) dst[name] = arr;
    }
}
} // namespace

void validate_hermiticity(TermSum const& terms, std::vector<SiteBasisPtr> const& sites,
                          FermionConvention convention, double tolerance,
                          std::size_t max_expanded_terms) {
    if (!std::isfinite(tolerance) || tolerance < 0. || !max_expanded_terms)
        throw std::invalid_argument("Hermiticity validation requires finite nonnegative tolerance and positive max_expanded_terms");
    if (convention != FermionConvention::strict && convention != FermionConvention::lattice)
        throw std::invalid_argument("unknown model fermion convention");
    auto full = validated_hamiltonian_terms(terms, sites, {}, convention);
    // Adjunction follows resolved same-site multiplication, rather than
    // reevaluating those products with a different floating-point association.
    auto dagger = validated_hamiltonian_terms(full.dag(), sites, {}, convention);
    if (exact_tensors(full) == exact_tensors(dagger)) return;

    // Homogeneous local parity makes each word fix its odd-site pattern.
    // Original and adjoint contributions share one invertible JW string,
    // which therefore factors out of the residual and magnitude comparisons.
    auto original = hermitian_words(full, max_expanded_terms);
    auto adjoint = hermitian_words(dagger, max_expanded_terms);
    Exact nonidentity_scale = 0;
    for (auto const* words : {&original, &adjoint})
        for (auto const& [word, coefficient] : *words)
            if (!word.empty()) nonidentity_scale = std::max(nonidentity_scale, coefficient.norm_squared());
    Exact const tolerance_squared = Exact(tolerance) * Exact(tolerance);
    auto check = [&](Word const& word, ExactComplex const& a, ExactComplex const& b) {
        Exact const scale = word.empty() ? nonidentity_scale
            : std::max(a.norm_squared(), b.norm_squared());
        if ((a - b).norm_squared() > tolerance_squared * scale)
            throw std::invalid_argument("Hamiltonian must be Hermitian");
    };
    for (auto const& [word, a] : original) {
        auto found = adjoint.find(word);
        check(word, a, found == adjoint.end() ? ExactComplex{} : found->second);
    }
    for (auto const& [word, b] : adjoint)
        if (!(original.count(word) != 0)) check(word, {}, b);
}

void validate_hermiticity(Hamiltonian const& model, double tolerance, std::size_t max_expanded_terms) {
    std::vector<SiteBasisPtr> sites;
    sites.reserve(static_cast<std::size_t>(model.num_sites()));
    for (int i = 0; i < model.num_sites(); ++i) sites.push_back(model.site_for(i));
    validate_hermiticity(model.expand(), sites, model.fermion_convention(), tolerance, max_expanded_terms);
}

TermSum compile_local_terms(TermSum const& terms, std::vector<SiteBasisPtr> const& sites,
                            std::vector<std::string> const& conserved_qns,
                            FermionConvention convention) {
    auto full = validated_hamiltonian_terms(terms, sites, conserved_qns, convention)
        .canonicalise(convention == FermionConvention::strict, 0.);
    if (!conserved_qns.empty()) {
        full = full.split_factors_by_qn_delta();
        std::vector<Term> kept;
        for (Term const& term : full.terms()) {
            Word transitions;
            // Every factor now has one exact charge difference. Recover it
            // from an actual entry; qn_deltas is floating-point metadata and
            // must not decide which side of a tolerance boundary gets kept.
            bool zero = false;
            for (auto const& factor : term.factors) {
                bool found = false;
                for (int row = 0; row < factor.op.dim && !found; ++row)
                    for (int column = 0; column < factor.op.dim && !found; ++column)
                        if (factor.op(row, column) != Complex{}) {
                            transitions.push_back({factor.site, row, column});
                            found = true;
                        }
                if (!found) { zero = true; break; }
            }
            bool drop = false;
            for (std::string const& q : conserved_qns) {
                if (!zero && changes_charge(transitions, sites, q)) { drop = true; break; }
            }
            if (!drop) kept.push_back(term);
        }
        full = TermSum(std::move(kept));
    }
    return full;
}

CompiledOperator compile_local_operator(TermSum const& terms,
    std::vector<SiteBasisPtr> const& sites, std::vector<std::string> const& conserved,
    FermionConvention convention) {
    if (convention != FermionConvention::strict && convention != FermionConvention::lattice)
        throw std::invalid_argument("unknown model fermion convention");
    if (sites.size() > std::size_t(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Hamiltonian site count exceeds the native index range");
    auto const compiled = compile_local_terms(terms, sites, conserved, convention);
    CompiledOperator result;
    result.site_dims.reserve(sites.size());
    for (auto const& site : sites) result.site_dims.push_back(site->dim);
    result.terms.reserve(compiled.terms().size());
    for (auto const& term : compiled.terms()) {
        auto& output = result.terms.emplace_back();
        output.coefficient = term.effective_coeff();
        output.jordan_wigner = jw_chain_sites(term.factors, static_cast<int>(sites.size()));
        output.factors.reserve(term.factors.size());
        for (auto const& factor : term.factors)
            output.factors.push_back({factor.site, factor.op.matrix});
    }
    return result;
}

CompiledOperator compile_hamiltonian(TermSum const& terms,
                                                std::vector<SiteBasisPtr> const& sites,
                                                std::vector<std::string> const& conserved_qns,
                                                FermionConvention convention) {
    if (convention != FermionConvention::strict && convention != FermionConvention::lattice)
        throw std::invalid_argument("unknown model fermion convention");
    if (sites.size() > std::size_t(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Hamiltonian site count exceeds the native index range");
    auto const canonical = validated_hamiltonian_terms(terms, sites, conserved_qns, convention)
        .canonicalise(convention == FermionConvention::strict, 0.);
    auto out = emit_operator(canonical, static_cast<int>(sites.size()),
                             convention == FermionConvention::strict);
    out.site_dims.reserve(sites.size());
    for (auto const& site : sites) out.site_dims.push_back(site->dim);
    return out;
}

CompiledOperator compile_operator(TermSum const& terms,
                                  std::vector<SiteBasisPtr> const& sites,
                                  FermionConvention convention) {
    if (convention != FermionConvention::strict && convention != FermionConvention::lattice)
        throw std::invalid_argument("unknown model fermion convention");
    if (sites.size() > std::size_t(std::numeric_limits<int>::max()))
        throw std::invalid_argument("operator site count exceeds the native index range");
    validate_local_inputs(terms, sites, {});
    bool const strict = convention == FermionConvention::strict;
    auto out = emit_operator(terms.canonicalise(strict, 0.), static_cast<int>(sites.size()), strict);
    out.site_dims.reserve(sites.size());
    for (auto const& site : sites) out.site_dims.push_back(site->dim);
    return out;
}

CompiledModel compile_hamiltonian(Hamiltonian const& model) {
    std::vector<SiteBasisPtr> sites;
    sites.reserve(std::size_t(model.num_sites()));
    for (int i = 0; i < model.num_sites(); ++i) sites.push_back(model.site_for(i));
    CompiledModel out{compile_hamiltonian(model.expand(), sites, model.conserved_qns(), model.fermion_convention())};
    add_model_metadata(out, model);
    return out;
}

CompiledModel Hamiltonian::to_local_decomposition() const {
    std::vector<SiteBasisPtr> sites;
    sites.reserve(std::size_t(num_sites()));
    for (int i = 0; i < num_sites(); ++i) sites.push_back(site_for(i));
    auto full = compile_local_terms(expand(), sites, conserved_qns_, convention_);
    return lower_operator(full, true);
}

CompiledModel Hamiltonian::operator_decomposition() const {
    return operator_decomposition(expand());
}

CompiledModel Hamiltonian::operator_decomposition(TermSum const& op) const {
    bool const strict = convention_ == FermionConvention::strict;
    return lower_operator(op.canonicalise(strict), strict);
}

CompiledModel Hamiltonian::lower_operator(TermSum const& full, bool jordan_wigner) const {
    CompiledModel out{emit_operator(full, num_sites(), jordan_wigner)};
    add_model_metadata(out, *this);
    return out;
}

} // namespace models
} // namespace alps
