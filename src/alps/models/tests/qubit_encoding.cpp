// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Qubit register layout of physical bases: widths, offsets and code words,
// bit and state-index round trips, statevector interchange and the rejected
// bases. Cases follow skilledwolf/alpsx (branch agent/v3-execution-and-api-audit
// at f8ee6e98b, tests/alps/test_spin_pauli.py and
// src/model/tests/models/pauli_export.cpp), adapted to GoogleTest on
// PhysicalBasis values.
#include <alps/models/builtin.hpp>
#include <alps/models/qubit_encoding.hpp>
#include <alps/models/site.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace lat = alps::lattice;
namespace mdl = alps::models;
using mdl::QubitEncoding;
using mdl::QubitEncodingScheme;
using Complex = std::complex<double>;

namespace {
template <class F>
void expect_invalid(F&& f, std::string const& fragment) {
    try {
        f();
        FAIL() << "expected invalid_argument containing: " << fragment;
    } catch (std::invalid_argument const& e) {
        EXPECT_NE(std::string(e.what()).find(fragment), std::string::npos) << e.what();
    }
}
} // namespace

TEST(QubitEncoding, WidthsOffsetsAndCodewords) {
    auto const basis = mdl::physical_basis({mdl::spin(.5), mdl::spin(1.), mdl::spin(2.), mdl::spin(7.5)});
    for (auto scheme : {QubitEncodingScheme::binary, QubitEncodingScheme::gray}) {
        QubitEncoding const encoding(basis, scheme);
        EXPECT_EQ(encoding.scheme(), scheme);
        EXPECT_EQ(encoding.num_qubits(), 10);
        EXPECT_EQ(encoding.site_dimensions(), (std::vector<int>{2, 3, 5, 16}));
        EXPECT_EQ(encoding.site_qubits(), (std::vector<int>{1, 2, 3, 4}));
        EXPECT_EQ(encoding.site_offsets(), (std::vector<int>{0, 1, 3, 6}));
        ASSERT_EQ(encoding.sites().size(), 4u);
        EXPECT_EQ(encoding.sites()[2].labels.size(), 5u);
        for (std::size_t site = 0; site < 4; ++site)
            for (int j = 0; j < encoding.site_dimensions()[site]; ++j) {
                auto const u = static_cast<std::size_t>(j);
                EXPECT_EQ(encoding.codeword(site, j), scheme == QubitEncodingScheme::binary ? u : u ^ (u >> 1));
            }
        EXPECT_THROW(encoding.codeword(2, 5), std::invalid_argument);
        EXPECT_THROW(encoding.codeword(4, 0), std::invalid_argument);
    }
    // Gray code words of the spin-2 site: 0, 1, 3, 2, 6.
    QubitEncoding const gray(basis, QubitEncodingScheme::gray);
    EXPECT_EQ(gray.codeword(2, 3), 2u);
    EXPECT_EQ(gray.codeword(2, 4), 6u);
}

TEST(QubitEncoding, BitsAndStateIndexRoundTrip) {
    auto const basis = mdl::physical_basis({mdl::spin(1.), mdl::spin(1.)});
    for (auto scheme : {QubitEncodingScheme::binary, QubitEncodingScheme::gray}) {
        QubitEncoding const encoding(basis, scheme);
        ASSERT_EQ(encoding.num_qubits(), 4);
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                std::vector<int> const physical{a, b};
                auto const bits = encoding.encode_basis_bits(physical);
                ASSERT_EQ(bits.size(), 4u);
                auto const decoded = encoding.decode_basis_bits(bits);
                ASSERT_TRUE(decoded.has_value());
                EXPECT_EQ(*decoded, physical);
                auto const code = [&](int j) { return scheme == QubitEncodingScheme::binary ? j : j ^ (j >> 1); };
                EXPECT_EQ(encoding.state_index(physical), static_cast<std::uint64_t>(code(a) + 4 * code(b)));
            }
        // The fourth code word of a three-state site is not physical.
        std::vector<int> const invalid = scheme == QubitEncodingScheme::binary ? std::vector<int>{1, 1, 0, 0}
                                                                               : std::vector<int>{0, 1, 0, 0};
        EXPECT_FALSE(encoding.decode_basis_bits(invalid).has_value());
        EXPECT_THROW(encoding.state_index({3, 0}), std::invalid_argument);
        EXPECT_THROW(encoding.state_index({0}), std::invalid_argument);
        EXPECT_THROW(encoding.encode_basis_bits({0, 0, 0}), std::invalid_argument);
        expect_invalid([&] { (void)encoding.decode_basis_bits({0, 0, 0, 2}); }, "binary bits");
        EXPECT_THROW(encoding.decode_basis_bits({0, 0, 0}), std::invalid_argument);
    }
}

TEST(QubitEncoding, StatevectorsProjectWithoutRenormalizing) {
    QubitEncoding const encoding(mdl::physical_basis({mdl::spin(1.), mdl::spin(1.)}));
    std::vector<Complex> psi(9);
    for (std::size_t k = 0; k < 9; ++k) psi[k] = Complex(.1 * static_cast<double>(k) - .3, .05 * static_cast<double>(k * k));
    auto const mapped = encoding.encode_statevector(psi);
    ASSERT_EQ(mapped.size(), 16u);
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
            EXPECT_EQ(mapped[encoding.state_index({a, b})], psi[static_cast<std::size_t>(3 * a + b)]);
    EXPECT_EQ(mapped[3], Complex{});
    EXPECT_EQ(mapped[12], Complex{});
    auto const decoded = encoding.decode_statevector(mapped);
    EXPECT_EQ(decoded.amplitudes, psi);
    EXPECT_EQ(decoded.discarded_norm_squared, 0.);
    EXPECT_NEAR(encoding.physical_weight(mapped), 1., 1e-15);

    std::vector<Complex> state(16);
    state[encoding.state_index({1, 0})] = Complex(0., .6);
    state[3] = .8; // code 3 on site 0: not physical
    auto const lossy = encoding.decode_statevector(state);
    EXPECT_NEAR(lossy.discarded_norm_squared, .64, 1e-15);
    double kept = 0.;
    for (auto z : lossy.amplitudes) kept += std::norm(z);
    EXPECT_NEAR(kept, .36, 1e-15);
    EXPECT_NEAR(encoding.physical_weight(state), .36, 1e-15);
    expect_invalid([&] { (void)encoding.physical_weight(std::vector<Complex>(16)); }, "positive norm");
    EXPECT_THROW(encoding.encode_statevector(psi, 15), std::length_error);
    EXPECT_THROW(encoding.encode_statevector(std::vector<Complex>(8)), std::invalid_argument);
    EXPECT_THROW(encoding.decode_statevector(std::vector<Complex>(9)), std::invalid_argument);
}

TEST(QubitEncoding, RejectsUnsupportedBases) {
    EXPECT_THROW(QubitEncoding(mdl::PhysicalBasis{}), std::invalid_argument);
    EXPECT_THROW(QubitEncoding(mdl::physical_basis({mdl::spin(.5)}), QubitEncodingScheme::binary, 0), std::invalid_argument);
    mdl::PhysicalBasis single;
    single.sites.push_back({{"only"}, mdl::Statistics::bosonic, {}, {1.}});
    expect_invalid([&] { QubitEncoding e(single); }, "dimension >= 2");
    expect_invalid([&] { QubitEncoding e(mdl::physical_basis({mdl::spinless_fermion()}), QubitEncodingScheme::gray); },
                   "binary dim-2/4");
    mdl::PhysicalBasis three;
    three.sites.push_back({{"0", "1", "2"}, mdl::Statistics::fermionic, {{"N", {0., 1., 2.}}}, {1., -1., 1.}});
    expect_invalid([&] { QubitEncoding e(three); }, "binary dim-2/4");
    mdl::PhysicalBasis parity;
    parity.sites.push_back({{"empty", "full"}, mdl::Statistics::fermionic, {{"N", {0., 1.}}}, {1., 1.}});
    expect_invalid([&] { QubitEncoding e(parity); }, "occupation parity");
    expect_invalid([&] { QubitEncoding e(mdl::physical_basis({mdl::spinless_fermion(), mdl::spin(1.)})); },
                   "separately validated");
    // Spin-1/2 and fermions mix; a 4x4 complex workspace needs 256 bytes.
    EXPECT_NO_THROW(QubitEncoding(mdl::physical_basis({mdl::spinless_fermion(), mdl::spin(.5)})));
    EXPECT_THROW(QubitEncoding(mdl::physical_basis({mdl::spin(1.)}), QubitEncodingScheme::binary, 255), std::length_error);
    EXPECT_NO_THROW(QubitEncoding(mdl::physical_basis({mdl::spin(1.)}), QubitEncodingScheme::binary, 256));
    // The basis itself is validated first.
    mdl::PhysicalBasis bad;
    bad.sites.push_back({{"up", "down"}, mdl::Statistics::bosonic, {{"Sz", {.5}}}, {1., 1.}});
    EXPECT_THROW(QubitEncoding e(bad), std::invalid_argument);
}

TEST(QubitEncoding, SnapshotMatchesTheModelBasis) {
    auto const model = mdl::heisenberg(lat::chain(2), 1., 0., 1.);
    auto const encoding = mdl::make_qubit_encoding(model);
    EXPECT_EQ(encoding.num_qubits(), 4);
    auto basis = mdl::physical_basis(model);
    EXPECT_NO_THROW(encoding.require_same_sites(basis));
    // Conserved charges and the convention are not part of the snapshot.
    basis.conserved_quantum_numbers.clear();
    basis.fermion_convention = mdl::FermionConvention::lattice;
    EXPECT_NO_THROW(encoding.require_same_sites(basis));
    auto reordered = basis;
    std::reverse(reordered.sites[0].quantum_numbers["Sz"].begin(), reordered.sites[0].quantum_numbers["Sz"].end());
    expect_invalid([&] { encoding.require_same_sites(reordered); }, "basis mismatch");
    auto relabelled = basis;
    std::reverse(relabelled.sites[1].labels.begin(), relabelled.sites[1].labels.end());
    expect_invalid([&] { encoding.require_same_sites(relabelled); }, "basis mismatch");
    auto longer = basis;
    longer.sites.push_back(basis.sites[0]);
    expect_invalid([&] { encoding.require_same_sites(longer); }, "count mismatch");
}
