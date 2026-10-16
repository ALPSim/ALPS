// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/model/src/operators/qubit_encoding.cpp. Adaptation: the encoding is built
// from and compared against a PhysicalBasis (validate_basis supplies the shape
// and finiteness checks alpsx repeated on SiteBasis), <alps/models/...>
// include names.
#include <alps/models/qubit_encoding.hpp>
#include <alps/models/hamiltonian.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace alps {
namespace models {
namespace {
using Complex = std::complex<double>;

double weight(Complex z) {
    double w = std::norm(z);
    if (!std::isfinite(w)) throw std::invalid_argument("state amplitudes must have finite squared norm");
    return w;
}
} // namespace

QubitEncoding::QubitEncoding(PhysicalBasis const& basis, QubitEncodingScheme scheme, std::size_t max_bytes)
    : scheme_(scheme) {
    if (basis.sites.empty() || max_bytes == 0)
        throw std::invalid_argument("encoding requires sites and positive max_local_matrix_bytes");
    if (scheme != QubitEncodingScheme::binary && scheme != QubitEncodingScheme::gray)
        throw std::invalid_argument("unknown qubit encoding scheme");
    validate_basis(basis);
    bool fermionic = false, higher_boson = false;
    for (auto const& s : basis.sites) {
        int const dim = static_cast<int>(s.labels.size());
        if (dim < 2) throw std::invalid_argument("encoding needs a local dimension >= 2");
        int width = 0;
        for (unsigned n = static_cast<unsigned>(dim - 1); n; n >>= 1) ++width;
        if (width >= std::numeric_limits<std::size_t>::digits)
            throw std::length_error("local encoding width exceeds native range");
        std::size_t const D = std::size_t{1} << width;
        if (D > max_bytes / sizeof(Complex) / D)
            throw std::length_error("encoded local matrix exceeds max_local_matrix_bytes");
        if (num_qubits_ > std::numeric_limits<int>::max() - width)
            throw std::length_error("qubit register width exceeds native range");
        bool const f = s.statistics == Statistics::fermionic;
        if (f && (scheme != QubitEncodingScheme::binary || (dim != 2 && dim != 4)))
            throw std::invalid_argument("fermionic encoding requires binary dim-2/4 occupation bases");
        for (int i = 0; i < dim; ++i) {
            double const expected = f && (i == 1 || i == 2) ? -1. : 1.;
            if (s.parity[static_cast<std::size_t>(i)] != expected)
                throw std::invalid_argument("encoding requires bosonic or binary occupation parity");
        }
        sites_.push_back(s);
        dimensions_.push_back(dim);
        offsets_.push_back(num_qubits_);
        widths_.push_back(width);
        num_qubits_ += width;
        fermionic |= f;
        higher_boson |= !f && dim > 2;
    }
    if (fermionic && higher_boson)
        throw std::invalid_argument("higher-dimensional spin/fermion mixtures require a separately validated encoding");
}

QubitEncoding make_qubit_encoding(Hamiltonian const& model, QubitEncodingScheme scheme, std::size_t bytes) {
    return QubitEncoding(physical_basis(model), scheme, bytes);
}

void QubitEncoding::require_same_sites(PhysicalBasis const& basis) const {
    if (basis.sites.size() != sites_.size()) throw std::invalid_argument("encoding basis/site count mismatch");
    for (std::size_t i = 0; i < sites_.size(); ++i)
        if (basis.sites[i] != sites_[i]) throw std::invalid_argument("encoding physical basis mismatch");
}

std::size_t QubitEncoding::codeword(std::size_t site, int local) const {
    if (site >= sites_.size() || local < 0 || local >= dimensions_[site])
        throw std::invalid_argument("invalid physical local basis index");
    auto const j = static_cast<std::size_t>(local);
    return scheme_ == QubitEncodingScheme::gray ? j ^ (j >> 1) : j;
}

std::vector<int> QubitEncoding::encode_basis_bits(std::vector<int> const& indices) const {
    if (indices.size() != sites_.size()) throw std::invalid_argument("expected one local index per site");
    std::vector<int> bits;
    bits.reserve(static_cast<std::size_t>(num_qubits_));
    for (std::size_t s = 0; s < sites_.size(); ++s) {
        auto const code = codeword(s, indices[s]);
        for (int b = 0; b < widths_[s]; ++b) bits.push_back(static_cast<int>((code >> b) & 1));
    }
    return bits;
}

std::optional<std::vector<int>> QubitEncoding::decode_basis_bits(std::vector<int> const& bits) const {
    if (bits.size() != static_cast<std::size_t>(num_qubits_)) throw std::invalid_argument("encoded bit count mismatch");
    for (int bit : bits) if (bit != 0 && bit != 1) throw std::invalid_argument("expected binary bits");
    std::vector<int> indices;
    for (std::size_t s = 0; s < sites_.size(); ++s) {
        std::size_t code = 0;
        for (int b = 0; b < widths_[s]; ++b)
            code |= static_cast<std::size_t>(bits[static_cast<std::size_t>(offsets_[s] + b)]) << b;
        std::size_t index = code;
        if (scheme_ == QubitEncodingScheme::gray)
            for (auto n = code >> 1; n; n >>= 1) index ^= n;
        if (index >= static_cast<std::size_t>(dimensions_[s])) return std::nullopt;
        indices.push_back(static_cast<int>(index));
    }
    return indices;
}

std::uint64_t QubitEncoding::state_index(std::vector<int> const& indices) const {
    if (num_qubits_ > 63) throw std::invalid_argument("state_index supports at most 63 qubits; use encode_basis_bits");
    if (indices.size() != sites_.size()) throw std::invalid_argument("expected one local index per site");
    std::uint64_t result = 0;
    for (std::size_t s = 0; s < sites_.size(); ++s)
        result |= static_cast<std::uint64_t>(codeword(s, indices[s])) << offsets_[s];
    return result;
}

std::pair<std::size_t, std::size_t> QubitEncoding::state_dimensions(std::size_t limit) const {
    if (!limit || num_qubits_ >= std::numeric_limits<std::size_t>::digits)
        throw std::length_error("state mapping exceeds max_amplitudes");
    std::size_t const encoded = std::size_t{1} << num_qubits_;
    if (encoded > limit) throw std::length_error("state mapping exceeds max_amplitudes");
    std::size_t physical = 1;
    for (int d : dimensions_) physical *= static_cast<std::size_t>(d); // bounded by encoded
    return {physical, encoded};
}

std::optional<std::size_t> QubitEncoding::physical_index(std::size_t encoded) const {
    std::size_t result = 0;
    for (std::size_t s = 0; s < sites_.size(); ++s) {
        auto const code = (encoded >> offsets_[s]) & ((std::size_t{1} << widths_[s]) - 1);
        auto index = code;
        if (scheme_ == QubitEncodingScheme::gray)
            for (auto n = code >> 1; n; n >>= 1) index ^= n;
        if (index >= static_cast<std::size_t>(dimensions_[s])) return std::nullopt;
        result = result * static_cast<std::size_t>(dimensions_[s]) + index;
    }
    return result;
}

std::vector<Complex> QubitEncoding::encode_statevector(std::vector<Complex> const& physical, std::size_t limit) const {
    auto const [d, D] = state_dimensions(limit);
    if (physical.size() != d) throw std::invalid_argument("physical state size mismatch");
    for (auto z : physical) (void)weight(z);
    std::vector<Complex> result(D);
    for (std::size_t i = 0; i < D; ++i) if (auto j = physical_index(i)) result[i] = physical[*j];
    return result;
}

DecodedQubitState QubitEncoding::decode_statevector(std::vector<Complex> const& encoded, std::size_t limit) const {
    auto const [d, D] = state_dimensions(limit);
    if (encoded.size() != D) throw std::invalid_argument("encoded state size mismatch");
    DecodedQubitState result;
    result.amplitudes.resize(d);
    for (std::size_t i = 0; i < D; ++i) {
        double const w = weight(encoded[i]);
        if (auto j = physical_index(i)) result.amplitudes[*j] = encoded[i];
        else result.discarded_norm_squared += w;
    }
    if (!std::isfinite(result.discarded_norm_squared)) throw std::invalid_argument("state norm overflow");
    return result;
}

double QubitEncoding::physical_weight(std::vector<Complex> const& encoded, std::size_t limit) const {
    auto const [d, D] = state_dimensions(limit);
    (void)d;
    if (encoded.size() != D) throw std::invalid_argument("encoded state size mismatch");
    double good = 0., total = 0.;
    for (std::size_t i = 0; i < D; ++i) {
        double const w = weight(encoded[i]);
        total += w;
        if (physical_index(i)) good += w;
    }
    if (!(total > 0.) || !std::isfinite(total)) throw std::invalid_argument("physical_weight requires finite positive norm");
    return good / total;
}

} // namespace models
} // namespace alps
