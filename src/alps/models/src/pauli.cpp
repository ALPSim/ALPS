// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/model/src/operators/pauli.cpp. Adaptation: the terms arrive as a
// CompiledOperator (canonical local factors with explicit Jordan-Wigner
// strings), so alpsx's per-factor basis checks, canonicalise_term and
// jw_chain_sites calls are replaced by a dimension check, in-order products of
// same-site factors and a definite-parity check on fermionic sites; local
// spaces, charges and the fermion convention come from the PhysicalBasis;
// <alps/models/...> include names.
#include <alps/models/pauli.hpp>
#include <alps/models/hamiltonian.hpp>
#include <alps/models/site.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>

namespace alps {
namespace models {
namespace {
using Expansion = std::map<std::string, Complex>;
constexpr char letters[] = "IXYZ";

bool finite(Complex z) { return std::isfinite(z.real()) && std::isfinite(z.imag()); }

void check_size(std::size_t n, std::size_t limit) {
    if (n > limit) throw std::length_error("Pauli expansion exceeds max_terms");
}

void validate_options(PauliExportOptions const& options) {
    double const validation = options.validation_tolerance.value_or(options.tolerance);
    if (!std::isfinite(options.tolerance) || options.tolerance < 0 || !options.max_terms
        || !options.max_local_matrix_bytes || !std::isfinite(validation) || validation < 0)
        throw std::invalid_argument("Pauli export requires finite nonnegative tolerance and positive resource limits");
}

std::size_t register_dimension(int width, std::size_t bytes) {
    if (width < 1 || width >= std::numeric_limits<std::size_t>::digits)
        throw std::length_error("local encoding width exceeds native range");
    auto const D = std::size_t{1} << width;
    if (D > bytes / sizeof(Complex) / D)
        throw std::length_error("encoded local matrix exceeds max_local_matrix_bytes");
    return D;
}

// In-place tensorized Pauli transform: each row/column bit pair maps to
// I,X,Y,Z. Each stage divides by two, giving Tr(P^dagger A)/2^width.
// O(width * 4^width) arithmetic, one local matrix, never a global matrix.
Expansion transform(Matrix matrix, int width, std::size_t limit) {
    auto const D = std::size_t{1} << width;
    for (int bit = 0; bit < width; ++bit) {
        auto const mask = std::size_t{1} << bit;
        for (std::size_t row = 0; row < D; ++row) if (!(row & mask))
            for (std::size_t col = 0; col < D; ++col) if (!(col & mask)) {
                auto const a = matrix[row * D + col] * .5;
                auto const b = matrix[row * D + (col | mask)] * .5;
                auto const c = matrix[(row | mask) * D + col] * .5;
                auto const d = matrix[(row | mask) * D + (col | mask)] * .5;
                matrix[row * D + col] = a + d;
                matrix[row * D + (col | mask)] = b + c;
                matrix[(row | mask) * D + col] = Complex(0., 1.) * (b - c);
                matrix[(row | mask) * D + (col | mask)] = a - d;
            }
    }
    Expansion result;
    for (std::size_t row = 0; row < D; ++row) for (std::size_t col = 0; col < D; ++col) {
        auto const coefficient = matrix[row * D + col];
        if (!finite(coefficient)) throw std::invalid_argument("non-finite Pauli coefficient");
        if (coefficient == Complex{}) continue; // no local tolerance pruning
        check_size(result.size() + 1, limit);
        std::string word(static_cast<std::size_t>(width), 'I');
        for (int bit = 0; bit < width; ++bit)
            word[static_cast<std::size_t>(width - 1 - bit)] = letters[2 * ((row >> bit) & 1) + ((col >> bit) & 1)];
        result.emplace(std::move(word), coefficient);
    }
    return result;
}

// Cache identical local matrices across sites within one compilation only.
// Bounded by both an estimated byte budget and 64 entries; keys retain exact
// finite matrix values. Eviction changes cost, never numerical meaning.
class LocalCompiler {
public:
    LocalCompiler(QubitEncoding const& encoding, PauliExportOptions const& options, PauliOperator& out)
        : encoding_(encoding), options_(options), out_(out) {}
    std::shared_ptr<Expansion const> operator()(Matrix const& matrix, std::size_t site) {
        int const dim = encoding_.site_dimensions()[site], width = encoding_.site_qubits()[site];
        auto const D = register_dimension(width, options_.max_local_matrix_bytes);
        auto const d = static_cast<std::size_t>(dim);
        if (matrix.size() != d * d) throw std::invalid_argument("local matrix shape mismatch");
        std::vector<double> key{double(dim), double(encoding_.scheme() == QubitEncodingScheme::gray)};
        key.reserve(2 + 2 * matrix.size());
        for (auto z : matrix) {
            if (!finite(z)) throw std::invalid_argument("non-finite local matrix after multiplication");
            key.push_back(z.real()); key.push_back(z.imag());
        }
        if (auto found = cache_.find(key); found != cache_.end()) {
            ++out_.local_cache_hits;
            return found->second;
        }
        Matrix encoded(D * D);
        // Linear, identity-preserving trace completion on unused states.
        if (D != d) {
            Complex mean = 0.;
            for (std::size_t i = 0; i < d; ++i) mean += matrix[i * d + i] / double(dim);
            for (std::size_t i = 0; i < D; ++i) encoded[i * D + i] = mean;
        }
        for (int r = 0; r < dim; ++r) for (int c = 0; c < dim; ++c)
            encoded[encoding_.codeword(site, r) * D + encoding_.codeword(site, c)] =
                matrix[static_cast<std::size_t>(r) * d + static_cast<std::size_t>(c)];
        auto result = std::make_shared<Expansion const>(transform(std::move(encoded), width, options_.max_terms));
        ++out_.local_decompositions;
        auto const cost = key.size() * sizeof(double)
                        + result->size() * (sizeof(Complex) + static_cast<std::size_t>(width) + 96);
        if (cost <= options_.max_local_matrix_bytes) {
            if (cache_.size() >= 64 || cache_bytes_ > options_.max_local_matrix_bytes - cost) {
                cache_.clear(); cache_bytes_ = 0;
            }
            cache_.emplace(std::move(key), result);
            cache_bytes_ += cost;
        }
        return result;
    }
private:
    QubitEncoding const& encoding_;
    PauliExportOptions const& options_;
    PauliOperator& out_;
    std::map<std::vector<double>, std::shared_ptr<Expansion const>> cache_;
    std::size_t cache_bytes_ = 0;
};

PauliOperator empty_operator(QubitEncoding const& encoding) {
    PauliOperator out;
    out.encoding = encoding;
    out.num_qubits = encoding.num_qubits();
    out.site_qubits = encoding.site_qubits();
    for (std::size_t s = 0; s < encoding.sites().size(); ++s)
        if (encoding.sites()[s].statistics == Statistics::fermionic)
            for (int b = out.site_qubits[s] - 1; b >= 0; --b)
                out.fermion_mode_qubits.push_back(encoding.site_offsets()[s] + b);
    return out;
}

void finish(PauliOperator& out, PauliExportOptions const& options) {
    double const validation = options.validation_tolerance.value_or(options.tolerance);
    for (auto it = out.terms.begin(); it != out.terms.end();) {
        auto& [word, coefficient] = *it;
        if (!finite(coefficient)) throw std::invalid_argument("non-finite accumulated Pauli coefficient");
        if (options.require_hermitian && std::abs(coefficient.imag()) > validation)
            throw std::invalid_argument("Hamiltonian must be Hermitian for Pauli export");
        int flips = 0;
        for (int q : out.fermion_mode_qubits) {
            char const p = word[static_cast<std::size_t>(out.num_qubits - 1 - q)];
            flips += p == 'X' || p == 'Y';
        }
        if (options.require_even && flips % 2 && std::abs(coefficient) > validation)
            throw std::invalid_argument("Hamiltonian must have even fermionic parity");
        if (std::abs(coefficient) <= options.tolerance) {
            out.discarded_l1 += std::abs(coefficient);
            it = out.terms.erase(it);
        } else {
            if (options.require_hermitian) {
                out.discarded_l1 += std::abs(coefficient.imag());
                coefficient = coefficient.real();
            }
            ++it;
        }
    }
    if (!std::isfinite(out.discarded_l1)) throw std::invalid_argument("discarded coefficient norm overflow");
}

std::pair<std::string, Complex> multiply(std::string a, std::string const& b) {
    Complex phase = 1.;
    for (std::size_t i = 0; i < a.size(); ++i) {
        char const x = a[i], y = b[i];
        if (y == 'I') continue;
        if (x == 'I') { a[i] = y; continue; }
        if (x == y) { a[i] = 'I'; continue; }
        bool const forward = (x == 'X' && y == 'Y') || (x == 'Y' && y == 'Z') || (x == 'Z' && y == 'X');
        phase *= Complex(0., forward ? 1. : -1.);
        a[i] = x != 'X' && y != 'X' ? 'X' : (x != 'Y' && y != 'Y' ? 'Y' : 'Z');
    }
    return {std::move(a), phase};
}

void validate_charges(PauliOperator const& out, PhysicalBasis const& basis, double tol, std::size_t limit,
                      LocalCompiler& compile) {
    for (auto const& name : basis.conserved_quantum_numbers) {
        Expansion charge;
        int offset = 0;
        for (std::size_t s = 0; s < basis.sites.size(); ++s) {
            auto const& site = basis.sites[s];
            auto const dim = site.labels.size();
            auto const it = site.quantum_numbers.find(name);
            if (it == site.quantum_numbers.end() || it->second.size() != dim)
                throw std::invalid_argument("missing site quantum number: " + name);
            Matrix matrix(dim * dim);
            for (std::size_t i = 0; i < dim; ++i) matrix[i * dim + i] = it->second[i];
            auto const parts = compile(matrix, s);
            for (auto const& [local, coeff] : *parts) {
                std::string word(static_cast<std::size_t>(out.num_qubits), 'I');
                word.replace(static_cast<std::size_t>(out.num_qubits - offset) - local.size(), local.size(), local);
                charge[word] += coeff;
                check_size(charge.size(), limit);
            }
            offset += out.site_qubits[s];
        }
        Expansion commutator;
        for (auto const& [h, hc] : out.terms) for (auto const& [q, qc] : charge) {
            auto [word, phase] = multiply(h, q);
            auto const reverse = multiply(q, h).second;
            if (phase != reverse) commutator[word] += hc * qc * (phase - reverse);
            check_size(commutator.size(), limit);
        }
        for (auto const& [word, coeff] : commutator) {
            (void)word;
            if (!finite(coeff) || std::abs(coeff) > tol)
                throw std::invalid_argument("Hamiltonian does not conserve declared quantum number " + name);
        }
    }
}

// Row-major product a*b of dim x dim matrices.
Matrix product(Matrix const& a, Matrix const& b, std::size_t dim) {
    Matrix out(dim * dim);
    for (std::size_t r = 0; r < dim; ++r)
        for (std::size_t k = 0; k < dim; ++k) {
            auto const x = a[r * dim + k];
            if (x == Complex{}) continue;
            for (std::size_t c = 0; c < dim; ++c) out[r * dim + c] += x * b[k * dim + c];
        }
    return out;
}

// Every nonzero entry of a fermionic local matrix must connect states of the
// same parity or states of opposite parity; a mixture has no CAR meaning.
void require_definite_parity(Matrix const& matrix, LocalBasis const& site) {
    auto const dim = site.labels.size();
    int sign = 0;
    for (std::size_t r = 0; r < dim; ++r) for (std::size_t c = 0; c < dim; ++c) {
        if (matrix[r * dim + c] == Complex{}) continue;
        int const s = site.parity[r] * site.parity[c] > 0 ? 1 : -1;
        if (sign && s != sign)
            throw std::invalid_argument("local matrix on a fermionic site mixes fermionic parities");
        sign = s;
    }
}
} // namespace

std::uint64_t PauliOperator::state_index(std::vector<int> const& indices) const {
    if (encoding) {
        if (num_qubits != encoding->num_qubits() || site_qubits != encoding->site_qubits())
            throw std::invalid_argument("inconsistent Pauli layout");
        return encoding->state_index(indices);
    }
    if (indices.size() != site_qubits.size() || num_qubits > 63)
        throw std::invalid_argument("state_index needs one local index per site and at most 63 qubits");
    std::uint64_t out = 0;
    int offset = 0;
    for (std::size_t s = 0; s < indices.size(); ++s) {
        int const width = site_qubits[s];
        if (width < 1 || width > 2 || offset + width > 63 || indices[s] < 0 || indices[s] >= (1 << width))
            throw std::invalid_argument("invalid encoded local basis index");
        out |= static_cast<std::uint64_t>(indices[s]) << offset;
        offset += width;
    }
    if (offset != num_qubits) throw std::invalid_argument("inconsistent Pauli layout");
    return out;
}

double PauliOperator::coefficient_l1(bool include_identity) const {
    double bound = 0.;
    for (auto const& [word, value] : terms)
        if (include_identity || word.find_first_not_of('I') != std::string::npos) bound += std::abs(value);
    if (!std::isfinite(bound)) throw std::invalid_argument("coefficient norm overflow");
    return bound;
}

PauliOperator export_pauli(Hamiltonian const& model, PauliExportOptions const& options) {
    return export_pauli(model, make_qubit_encoding(model, QubitEncodingScheme::binary, options.max_local_matrix_bytes),
                        options);
}

PauliOperator export_pauli(Hamiltonian const& model, QubitEncoding const& encoding, PauliExportOptions const& options) {
    return export_pauli(compile_hamiltonian(model).op, physical_basis(model), encoding, options);
}

PauliOperator export_pauli(CompiledOperator const& op, PhysicalBasis const& basis, PauliExportOptions const& options) {
    return export_pauli(op, basis, QubitEncoding(basis, QubitEncodingScheme::binary, options.max_local_matrix_bytes),
                        options);
}

PauliOperator export_pauli(CompiledOperator const& op, PhysicalBasis const& basis, QubitEncoding const& encoding,
                           PauliExportOptions const& options) {
    validate_options(options);
    validate_basis(basis);
    encoding.require_same_sites(basis);
    if (op.site_dims.size() != basis.sites.size())
        throw std::invalid_argument("compiled operator site count differs from the physical basis");
    for (std::size_t s = 0; s < basis.sites.size(); ++s)
        if (op.site_dims[s] < 0 || static_cast<std::size_t>(op.site_dims[s]) != basis.sites[s].labels.size())
            throw std::invalid_argument("compiled operator local dimension differs from the physical basis at site "
                                        + std::to_string(s));
    for (int width : encoding.site_qubits()) (void)register_dimension(width, options.max_local_matrix_bytes);
    auto out = empty_operator(encoding);
    auto const& offsets = encoding.site_offsets();
    if (!out.fermion_mode_qubits.empty() && basis.fermion_convention != FermionConvention::strict)
        throw std::invalid_argument("fermionic Pauli export requires strict convention");
    LocalCompiler compile(encoding, options, out);

    for (auto const& term : op.terms) {
        if (!finite(term.coefficient)) throw std::invalid_argument("non-finite coefficient");
        std::map<int, Matrix> local;
        for (auto const& factor : term.factors) {
            if (factor.site < 0 || static_cast<std::size_t>(factor.site) >= basis.sites.size())
                throw std::invalid_argument("Pauli factor site is out of range");
            auto const& site = basis.sites[static_cast<std::size_t>(factor.site)];
            auto const dim = site.labels.size();
            if (factor.matrix.size() != dim * dim) throw std::invalid_argument("Pauli factor dimension mismatch");
            for (auto value : factor.matrix)
                if (!finite(value)) throw std::invalid_argument("non-finite local matrix");
            if (site.statistics == Statistics::fermionic) require_definite_parity(factor.matrix, site);
            auto [it, inserted] = local.try_emplace(factor.site, factor.matrix);
            if (!inserted) it->second = product(it->second, factor.matrix, dim);
        }
        for (int s : term.jordan_wigner) {
            if (s < 0 || static_cast<std::size_t>(s) >= basis.sites.size())
                throw std::invalid_argument("Jordan-Wigner site is out of range");
            auto const& site = basis.sites[static_cast<std::size_t>(s)];
            auto const dim = site.labels.size();
            auto [it, inserted] = local.try_emplace(s, dim * dim, Complex{});
            if (inserted) for (std::size_t i = 0; i < dim; ++i) it->second[i * dim + i] = 1.;
            for (std::size_t r = 0; r < dim; ++r) for (std::size_t c = 0; c < dim; ++c)
                it->second[r * dim + c] *= site.parity[c];
        }
        Expansion expansion{{std::string(static_cast<std::size_t>(out.num_qubits), 'I'), term.coefficient}};
        for (auto const& [s, matrix] : local) {
            auto const parts_owner = compile(matrix, static_cast<std::size_t>(s));
            auto const& parts = *parts_owner;
            if (!parts.empty() && expansion.size() > options.max_terms / parts.size())
                throw std::length_error("Pauli expansion exceeds max_terms");
            Expansion next;
            for (auto const& [word, coefficient] : expansion) for (auto const& [part, weight] : parts) {
                auto key = word;
                key.replace(static_cast<std::size_t>(out.num_qubits - offsets[static_cast<std::size_t>(s)]) - part.size(),
                            part.size(), part);
                next.emplace(std::move(key), coefficient * weight);
            }
            expansion = std::move(next);
        }
        for (auto const& [word, coefficient] : expansion) {
            out.terms[word] += coefficient;
            if (!finite(out.terms[word])) throw std::invalid_argument("non-finite accumulated Pauli coefficient");
            check_size(out.terms.size(), options.max_terms);
        }
    }
    // Validate the unpruned operator, so pruning cannot conceal a stale sector.
    if (options.require_conservation)
        validate_charges(out, basis, options.validation_tolerance.value_or(options.tolerance), options.max_terms, compile);
    finish(out, options);
    return out;
}

PauliOperator local_valid_projector(QubitEncoding const& encoding, std::size_t site, PauliExportOptions const& options) {
    validate_options(options);
    if (site >= encoding.site_dimensions().size()) throw std::invalid_argument("projector site out of range");
    int const width = encoding.site_qubits()[site];
    auto const D = register_dimension(width, options.max_local_matrix_bytes);
    Matrix projector(D * D);
    for (int i = 0; i < encoding.site_dimensions()[site]; ++i) {
        auto const code = encoding.codeword(site, i);
        projector[code * D + code] = 1.;
    }
    auto out = empty_operator(encoding);
    for (auto const& [local, value] : transform(std::move(projector), width, options.max_terms)) {
        std::string word(static_cast<std::size_t>(out.num_qubits), 'I');
        word.replace(static_cast<std::size_t>(out.num_qubits - encoding.site_offsets()[site]) - local.size(),
                     local.size(), local);
        out.terms.emplace(std::move(word), value);
    }
    finish(out, options);
    return out;
}

PauliOperator invalid_state_penalty(QubitEncoding const& encoding, double strength, PauliExportOptions const& options) {
    validate_options(options);
    if (!std::isfinite(strength) || strength < 0.)
        throw std::invalid_argument("penalty strength must be finite and nonnegative");
    auto out = empty_operator(encoding);
    // Prune only after scaling and summing; otherwise the reported bound
    // would not cover a projector's discarded terms amplified by strength.
    auto exact = options;
    exact.tolerance = 0.;
    for (std::size_t s = 0; s < encoding.site_dimensions().size(); ++s) {
        out.terms[std::string(static_cast<std::size_t>(out.num_qubits), 'I')] += strength;
        for (auto const& [word, value] : local_valid_projector(encoding, s, exact).terms)
            out.terms[word] -= strength * value;
        check_size(out.terms.size(), options.max_terms);
    }
    finish(out, options);
    return out;
}

} // namespace models
} // namespace alps
