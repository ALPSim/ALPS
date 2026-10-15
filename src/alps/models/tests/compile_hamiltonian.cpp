// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/tests/models/compile_hamiltonian.cpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports;
// dense assembly replaces the ED kernel of the original.
#include <alps/models/builtin.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstdint>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
namespace mdl = alps::models;
using C = mdl::Complex;
using Input = alps::models::CompiledOperator;

// Consume the public compiled descriptor by dense Kronecker assembly: per
// term, the product over all sites of (factor op [times F when the site is on
// the term's JW string] | bare F on factor-less JW sites | identity), with site
// zero most significant. Expected matrix entries below are specified
// independently from spin algebra or ordered Fock-state signs.
std::vector<C> matrix(Input const& input, std::vector<std::vector<double>> const& parity = {}) {
    std::size_t const n = input.site_dims.size();
    std::size_t dimension = 1;
    for (int dim : input.site_dims) dimension *= static_cast<std::size_t>(dim);
    std::vector<C> result(dimension * dimension);
    for (auto const& term : input.terms) {
        std::vector<C> acc{1.};
        std::size_t acc_dim = 1;
        for (std::size_t s = 0; s < n; ++s) {
            std::size_t const sd = static_cast<std::size_t>(input.site_dims[s]);
            std::vector<C> local(sd * sd);
            bool const in_jw = std::find(term.jordan_wigner.begin(), term.jordan_wigner.end(),
                                         static_cast<int>(s)) != term.jordan_wigner.end();
            auto const fit = std::find_if(term.factors.begin(), term.factors.end(),
                [s](auto const& factor) { return factor.site == static_cast<int>(s); });
            auto const f = [&](std::size_t c) { return in_jw ? (parity.empty() ? 1. : parity[s][c]) : 1.; };
            if (fit != term.factors.end()) {
                for (std::size_t r = 0; r < sd; ++r)
                    for (std::size_t c = 0; c < sd; ++c) local[r * sd + c] = fit->matrix[r * sd + c] * f(c);
            } else {
                for (std::size_t r = 0; r < sd; ++r) local[r * sd + r] = f(r);
            }
            std::vector<C> next(acc_dim * sd * acc_dim * sd);
            for (std::size_t r = 0; r < acc_dim; ++r)
                for (std::size_t c = 0; c < acc_dim; ++c) {
                    if (acc[r * acc_dim + c] == 0.) continue;
                    for (std::size_t i = 0; i < sd; ++i)
                        for (std::size_t j = 0; j < sd; ++j)
                            next[(r * sd + i) * acc_dim * sd + c * sd + j] = acc[r * acc_dim + c] * local[i * sd + j];
                }
            acc = std::move(next);
            acc_dim *= sd;
        }
        for (std::size_t i = 0; i < result.size(); ++i) result[i] += term.coefficient * acc[i];
    }
    return result;
}
} // namespace

TEST(CompiledOperator, OwnsChargeChangingObservableWithoutModelMetadata) {
    mdl::CompiledOperator observable;
    C const amplitude(.75, .25);
    {
        auto writable = std::make_shared<mdl::SiteBasis>(*mdl::spin());
        auto terms = amplitude * mdl::TermSum::located(writable, "Sp", 0);
        observable = mdl::compile_operator(terms, {writable});
        writable->operators["Sp"].matrix.assign(4, 19.);
    }
    // This non-Hermitian, charge-changing observable remains usable after
    // all source algebra and local basis objects have been destroyed.
    EXPECT_EQ(matrix(observable), (std::vector<C>{0., amplitude, 0., 0.}));
    auto independent = observable;
    ASSERT_EQ(independent.terms.size(), 1);
    ASSERT_EQ(independent.terms[0].factors.size(), 1);
    independent.terms[0].factors[0].matrix[1] = 7.;
    EXPECT_EQ(matrix(observable), (std::vector<C>{0., amplitude, 0., 0.}));
}

TEST(NativeObservableLowering, ExternalScalePreservesNearIdentityCorrections) {
    for (bool diagonal : {false, true}) {
        mdl::SiteBasis description = *mdl::spin();
        double const scale = diagonal ? 0x1p30 : 0x1p45;
        mdl::Matrix const local = diagonal ? mdl::Matrix{1., 0., 0., 1. + 1. / scale}
                                           : mdl::Matrix{1., 1. / scale, 1. / scale, 1.};
        description.operators["A"] = {"A", 2, local};
        auto site = mdl::finalize_site(std::move(description));
        auto const compiled = mdl::compile_operator(
            scale * mdl::TermSum::located(site, "A", 0), {site});
        std::vector<C> const expected = diagonal ? std::vector<C>{scale, 0., 0., scale + 1.}
                                                : std::vector<C>{scale, 1., 1., scale};
        EXPECT_EQ(matrix(compiled), expected);
    }
}

TEST(NativeObservableLowering, CancellingCoefficientsRetainTheUnitRemainder) {
    auto site = mdl::spin();
    auto const sz = mdl::TermSum::located(site, "Sz", 0);
    for (auto const& coefficients : std::vector<std::vector<C>>{
             {0x1p54, 1., -0x1p54}, {0x1p54, -0x1p54, 1.}}) {
        mdl::TermSum observable;
        for (auto coefficient : coefficients) observable = observable + coefficient * sz;
        auto const compiled = mdl::compile_operator(observable, {site});
        EXPECT_EQ(matrix(compiled), (std::vector<C>{.5, 0., 0., -.5}));
    }
}

TEST(NativeObservableLowering, OddFermionCreationKeepsTheJordanWignerSigns) {
    auto site = mdl::spinless_fermion();
    C const amplitude(.75, .25);
    auto const compiled = mdl::compile_operator(
        amplitude * mdl::TermSum::located(site, "c_dag", 2), {site, site, site});
    std::vector<C> expected(64);
    // Creation at site 2 crosses the two earlier occupation bits. Product
    // rows are |000>, |001>, ..., |111>, with site zero most significant.
    expected[1 * 8 + 0] = expected[7 * 8 + 6] = amplitude;
    expected[3 * 8 + 2] = expected[5 * 8 + 4] = -amplitude;
    EXPECT_EQ(matrix(compiled, {{1., -1.}, {1., -1.}, {1., -1.}}), expected);
}

TEST(NativeObservableLowering, NonHermitianSameSiteProductsKeepTheirOrder) {
    auto site = mdl::spin();
    auto const raising = mdl::TermSum::located(site, "Sp", 0);
    auto const sz = mdl::TermSum::located(site, "Sz", 0);
    EXPECT_EQ(matrix(mdl::compile_operator(raising * sz, {site})),
              (std::vector<C>{0., -.5, 0., 0.}));
    EXPECT_EQ(matrix(mdl::compile_operator(sz * raising, {site})),
              (std::vector<C>{0., .5, 0., 0.}));
}

TEST(CompileHamiltonian, RejectsFalseGlobalConservationWithoutChangingArbitraryOperators) {
    auto site = mdl::spin();
    auto op = [&](char const* name, int i) { return mdl::TermSum::located(site, name, i); };
    for (bool invariant_zero_sector : {false, true}) {
        mdl::Hamiltonian h(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Sz"});
        auto terms = invariant_zero_sector ? op("Sx", 0) * op("Sx", 1) - op("Sy", 0) * op("Sy", 1)
                                          : op("Sx", 0);
        h.add(terms);
        EXPECT_THROW(mdl::compile_hamiltonian(h), std::invalid_argument);
        EXPECT_THROW(mdl::compile_hamiltonian(terms, {site, site}, {"Sz"}), std::invalid_argument);
        // The pair-flip model changes Sz only outside the Sz=0 sector. Global
        // model validation still rejects it; intentional operator lowering
        // retains its all-up/all-down amplitude instead of projecting it.
        auto arbitrary = matrix(h.operator_decomposition().op);
        auto const target = invariant_zero_sector ? 3 : 2;
        EXPECT_EQ(arbitrary[target], C(.5));
        EXPECT_EQ(arbitrary[4 * target], C(.5));
    }
}

TEST(CompileHamiltonian, KeepsUnsplitCartesianTermsAndTheirComplexMatrices) {
    auto site = mdl::spin();
    auto op = [&](char const* name, int i) { return mdl::TermSum::located(site, name, i); };
    auto terms = op("Sx", 0) * op("Sx", 1) + op("Sy", 0) * op("Sy", 1);
    auto raw = mdl::compile_hamiltonian(terms, {site, site}, {"Sz"});
    ASSERT_EQ(raw.terms.size(), 2);
    for (std::size_t term = 0; term < 2; ++term) {
        ASSERT_EQ(raw.terms[term].factors.size(), 2);
        for (int index = 0; index < 2; ++index) {
            EXPECT_EQ(raw.terms[term].factors[index].site, index);
            EXPECT_EQ(raw.terms[term].factors[index].matrix,
                      site->op(term == 0 ? "Sx" : "Sy").matrix);
        }
    }
    std::vector<C> expected(16);
    expected[1 * 4 + 2] = expected[2 * 4 + 1] = .5;
    EXPECT_EQ(matrix(raw), expected);

    mdl::Hamiltonian h(alps::lattice::chain(2, 1., alps::lattice::Boundary::open), site,
                       mdl::FermionConvention::strict, {"Sz"});
    h.add(terms);
    auto compiled = mdl::compile_hamiltonian(h);
    EXPECT_EQ(matrix(compiled.op), expected);
    EXPECT_EQ(compiled.site_types, h.lat().site_types());
    ASSERT_EQ(compiled.edges.size(), 1);
    EXPECT_EQ(compiled.edges[0].src, 0);
    EXPECT_EQ(compiled.edges[0].dst, 1);
    EXPECT_EQ(compiled.diag_observables.at(0).at("Sz"), site->quantum_numbers.at("Sz"));
}

TEST(CompileHamiltonian, FullOperatorCrossSupportCancellationPrecedesAnyProjection) {
    mdl::SiteBasis basis = *mdl::spin();
    basis.operators["Pup"] = {"Pup", 2, {1., 0., 0., 0.}};
    basis.operators["Pdown"] = {"Pdown", 2, {0., 0., 0., 1.}};
    auto site = mdl::finalize_site(std::move(basis));
    auto op = [&](char const* name, int i) { return mdl::TermSum::located(site, name, i); };
    for (bool reverse : {false, true})
        for (double remainder : {0., 1.}) {
            std::vector<mdl::TermSum> pieces{0x1p50 * op("Sx", 0) * op("Pup", 1),
                0x1p50 * op("Sx", 0) * op("Pdown", 1), -(0x1p50 - remainder) * op("Sx", 0)};
            if (reverse) std::reverse(pieces.begin(), pieces.end());
            mdl::Hamiltonian h(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Sz"});
            for (auto const& piece : pieces) h.add(piece);
            h.add(op("Sz", 0) * op("Sz", 1));
            if (remainder != 0.) {
                EXPECT_THROW(mdl::compile_hamiltonian(h), std::invalid_argument);
            } else {
                auto compiled = mdl::compile_hamiltonian(h);
                EXPECT_EQ(compiled.op.terms.size(), 4); // Preserve unsplit contributions.
                std::vector<C> expected(16);
                expected[0] = expected[15] = .25;
                expected[5] = expected[10] = -.25;
                EXPECT_EQ(matrix(compiled.op), expected);
            }
        }
}

TEST(CompileHamiltonian, ValidationPrecedesAccurateLikeTermAccumulation) {
    auto site = mdl::spin();
    auto op = [&](char const* name) { return mdl::TermSum::located(site, name, 0); };
    for (auto const& coefficients : std::vector<std::vector<C>>{
             {0x1p54, 1., -0x1p54}, {0x1p54, -0x1p54, 1.}}) {
        mdl::TermSum broken, conserving;
        for (auto coefficient : coefficients) {
            broken = broken + coefficient * op("Sx");
            conserving = conserving + coefficient * op("Sz");
        }
        EXPECT_THROW(mdl::compile_hamiltonian(broken, {site}, {"Sz"}), std::invalid_argument);
        // The final exactly representable remainder is independent of order.
        auto const compiled = mdl::compile_hamiltonian(conserving, {site}, {"Sz"});
        ASSERT_EQ(compiled.terms.size(), 1);
        EXPECT_EQ(compiled.terms[0].coefficient, C(1.));
    }
}

TEST(CompileHamiltonian, NearIdentityEntriesArePreservedBeforeExternalScaling) {
    mdl::SiteBasis basis = *mdl::spin();
    basis.operators["NearI"] = {"NearI", 2, {1., 0x1p-45, 0x1p-45, 1.}};
    auto site = mdl::finalize_site(std::move(basis));
    auto terms = 0x1p45 * mdl::TermSum::located(site, "NearI", 0);
    auto compiled = mdl::compile_hamiltonian(terms, {site});
    EXPECT_EQ(matrix(compiled), (std::vector<C>{0x1p45, 1., 1., 0x1p45}));
    EXPECT_THROW(mdl::compile_hamiltonian(terms, {site}, {"Sz"}), std::invalid_argument);
}

TEST(CompileHamiltonian, ComplexHoppingRetainsBothFermionConventions) {
    auto site = mdl::spinless_fermion();
    C const amplitude(.75, .25);
    auto terms = amplitude * mdl::TermSum::located(site, "c_dag", 2)
                           * mdl::TermSum::located(site, "c", 0);
    terms = terms.add_hc();
    for (auto convention : {mdl::FermionConvention::strict, mdl::FermionConvention::lattice}) {
        auto compiled = mdl::compile_hamiltonian(terms, {site, site, site}, {"N"}, convention);
        auto actual = matrix(compiled, std::vector<std::vector<double>>(3, site->parity_operator));
        std::vector<C> expected(64);
        // Ordered Fock states: hopping 0->2 crosses the occupied middle
        // site only for |110> -> |011>, producing a minus sign in strict
        // convention. The lattice tensor convention has no such sign.
        C const crossed = convention == mdl::FermionConvention::strict ? -amplitude : amplitude;
        expected[1 * 8 + 4] = amplitude;
        expected[4 * 8 + 1] = std::conj(amplitude);
        expected[3 * 8 + 6] = crossed;
        expected[6 * 8 + 3] = std::conj(crossed);
        EXPECT_EQ(actual, expected);
    }
}

TEST(CompileHamiltonian, NeutralHeterogeneousChargeRowsRetainKondoExchange) {
    alps::lattice::UnitCell cell(1);
    cell.add_vertex({0.}, 0);
    cell.add_vertex({.5}, 1);
    cell.add_edge(0, 1, 1);
    auto h = mdl::kondo_lattice(cell.tile({1}), 0., 1.);
    auto compiled = mdl::compile_hamiltonian(h);
    EXPECT_EQ(compiled.op.site_dims, (std::vector<int>{4, 2}));
    EXPECT_EQ(h.site_for(1)->quantum_numbers.count("N"), 0);
    auto actual = matrix(compiled.op, {h.site_for(0)->parity_operator, h.site_for(1)->parity_operator});
    // Electron local rows are empty, up, down, doubly occupied. The empty
    // and double states carry no spin; the singly occupied 4x4 block is
    // Sz.s_z + (S+.s_- + S-.s_+)/2.
    std::vector<C> expected(64);
    expected[2 * 8 + 2] = expected[5 * 8 + 5] = .25;
    expected[3 * 8 + 3] = expected[4 * 8 + 4] = -.25;
    expected[3 * 8 + 4] = expected[4 * 8 + 3] = .5;
    EXPECT_EQ(actual, expected);
}

TEST(CompileHamiltonian, DescriptorOwnsMatricesAndModelMetadata) {
    auto site = std::make_shared<mdl::SiteBasis>(*mdl::spin());
    mdl::Hamiltonian h(alps::lattice::chain(2), site, mdl::FermionConvention::strict, {"Sz"});
    h.add(mdl::TermSum::located(site, "Sz", 0));
    auto compiled = mdl::compile_hamiltonian(h);
    auto before = matrix(compiled.op);
    site->operators["Sz"].matrix[0] = 19.;
    site->quantum_numbers["Sz"][0] = 19.;
    EXPECT_EQ(matrix(compiled.op), before);
    EXPECT_EQ(compiled.diag_observables.at(0).at("Sz"), (std::vector<double>{.5, -.5}));
}

TEST(CompileHamiltonian, ValidationDoesNotEnumerateTheManyBodyHilbertSpace) {
    auto h = mdl::heisenberg(alps::lattice::chain(64));
    auto compiled = mdl::compile_hamiltonian(h);
    EXPECT_EQ(compiled.op.site_dims, std::vector<int>(64, 2));
    EXPECT_FALSE(compiled.op.terms.empty());
    for (auto const& term : compiled.op.terms) EXPECT_LE(term.factors.size(), 2);
}

TEST(CompileHamiltonian, RejectsMalformedNumericalInputs) {
    auto site = mdl::spin();
    auto terms = mdl::TermSum::located(site, "Sz", 0);
    EXPECT_THROW(mdl::compile_hamiltonian(terms, {site}, {"unknown"}), std::invalid_argument);
    EXPECT_THROW(mdl::compile_hamiltonian(terms, {nullptr}, {"Sz"}), std::invalid_argument);
    EXPECT_THROW(mdl::compile_hamiltonian(terms, {site}, {}, static_cast<mdl::FermionConvention>(99)),
                 std::invalid_argument);
    EXPECT_THROW(mdl::compile_hamiltonian(std::numeric_limits<double>::infinity() * terms, {site}),
                 std::invalid_argument);
}
