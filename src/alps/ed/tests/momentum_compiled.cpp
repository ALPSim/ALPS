// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/tests/momentum_compiled.cpp.
// Adaptation: C++17 (no std::popcount, std::span or std::numbers),
// <alps/ed/...> include names; the dispatch cases follow the IETL driver's
// policy (Lanczos above dimension 1000 and 8*count + 200, SM always dense;
// defaults tolerance 1e-10 and min(10*dimension, 1000) steps).
#include <alps/ed/eigensolver.hpp>
#include <alps/ed/projected_basis.hpp>
#include <alps/ed/spectrum.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using Complex = std::complex<double>;
double const pi = std::acos(-1.);
int popcount(unsigned bits) {
    int count = 0;
    for (; bits; bits &= bits - 1) ++count;
    return count;
}
struct Problem {
    std::vector<std::int64_t> offsets{0}, columns, states;
    std::vector<Complex> values;
    std::vector<alps::ed::SpectrumLocalSpace> spaces;
    std::vector<alps::ed::SpectrumTranslation> translations;
    std::vector<alps::ed::SpectrumBlock> solve(int k = 128) const {
        return alps::ed::symmetry_spectrum(offsets, columns, values, offsets.size() - 1, states, spaces,
                                           translations, k);
    }
};

// Independent Fock-space reference: apply c_i^dagger c_j with the signs
// from annihilating and then creating in the ordered occupied orbitals.
Problem ring(int length, int particles, double phase) {
    Problem p;
    p.spaces.resize(static_cast<std::size_t>(length), {{1., -1.}, {{0., 1.}}});
    std::vector<std::size_t> translation(static_cast<std::size_t>(length));
    for (int i = 0; i < length; ++i) translation[static_cast<std::size_t>(i)] = static_cast<std::size_t>((i + 1) % length);
    p.translations.push_back({translation, static_cast<std::size_t>(length)});
    std::vector<unsigned> basis;
    for (unsigned bits = 0; bits < (1u << length); ++bits)
        if (popcount(bits) == particles) {
            basis.push_back(bits);
            for (int i = 0; i < length; ++i) p.states.push_back((bits >> i) & 1);
        }
    std::vector<Complex> matrix(basis.size() * basis.size());
    for (std::size_t col = 0; col < basis.size(); ++col)
        for (int i = 0; i < length; ++i)
            for (int reverse = 0; reverse < 2; ++reverse) {
                int const j = (i + 1) % length;
                int const source = reverse ? i : j, target = reverse ? j : i;
                unsigned bits = basis[col];
                if (!(bits & (1u << source)) || (bits & (1u << target))) continue;
                int sign = popcount(bits & ((1u << source) - 1));
                bits ^= 1u << source;
                sign += popcount(bits & ((1u << target) - 1));
                bits |= 1u << target;
                auto const row = static_cast<std::size_t>(std::find(basis.begin(), basis.end(), bits) - basis.begin());
                matrix[row * basis.size() + col] += (sign % 2 ? 1. : -1.) * std::polar(1., reverse ? -phase : phase);
            }
    for (std::size_t row = 0; row < basis.size(); ++row) {
        for (std::size_t col = 0; col < basis.size(); ++col)
            if (matrix[row * basis.size() + col] != Complex{}) {
                p.columns.push_back(static_cast<std::int64_t>(col));
                p.values.push_back(matrix[row * basis.size() + col]);
            }
        p.offsets.push_back(static_cast<std::int64_t>(p.values.size()));
    }
    return p;
}
// One uncharged block of the given dimension over binary site digits: a
// tridiagonal matrix with the given diagonal and hopping.
Problem tridiagonal(int dimension, Complex hopping, double slope, double offset) {
    Problem p;
    int sites = 1;
    while ((1 << sites) < dimension) ++sites;
    p.spaces.resize(static_cast<std::size_t>(sites), {{1., 1.}, {}});
    std::vector<std::size_t> identity(static_cast<std::size_t>(sites));
    std::iota(identity.begin(), identity.end(), std::size_t{0});
    p.translations.push_back({identity, 1});
    for (int row = 0; row < dimension; ++row) {
        for (int site = 0; site < sites; ++site) p.states.push_back((row >> site) & 1);
        if (row > 0) { p.columns.push_back(row - 1); p.values.push_back(hopping); }
        p.columns.push_back(row);
        p.values.emplace_back(offset + slope * row);
        if (row + 1 < dimension) { p.columns.push_back(row + 1); p.values.push_back(std::conj(hopping)); }
        p.offsets.push_back(static_cast<std::int64_t>(p.values.size()));
    }
    return p;
}
} // namespace

TEST(CompiledMomentum, ComplexHoppingRetainsTheOrientedDispersion) {
    auto p = ring(5, 1, .3);
    auto result = p.solve();
    ASSERT_EQ(result.size(), 5);
    for (auto const& block : result) {
        ASSERT_EQ(block.eigenvalues.size(), 1);
        EXPECT_DOUBLE_EQ(block.quantum_numbers[0], 1.);
        EXPECT_NEAR(block.eigenvalues[0], -2. * std::cos(block.momentum[0] - .3), 1e-13);
    }
}

TEST(CompiledMomentum, FermionicOrbitsMatchSumsOfDistinctOneParticleLevels) {
    int const length = 6;
    double const phase = .31;
    auto const result = ring(length, 2, phase).solve();
    ASSERT_EQ(result.size(), length);
    for (int k = 0; k < length; ++k) {
        std::vector<double> expected;
        for (int a = 0; a < length; ++a)
            for (int b = a + 1; b < length; ++b)
                if ((a + b) % length == k)
                    expected.push_back(-2. * std::cos(2. * pi * a / length - phase) - 2. * std::cos(2. * pi * b / length - phase));
        std::sort(expected.begin(), expected.end());
        auto const& block = result[static_cast<std::size_t>(k)];
        ASSERT_EQ(block.eigenvalues.size(), expected.size());
        for (std::size_t i = 0; i < expected.size(); ++i) EXPECT_NEAR(block.eigenvalues[i], expected[i], 1e-12);
    }
}

TEST(CompiledMomentum, BothSparseScalarTypesUseTheSharedBlockSolver) {
    // A 64-site chain is solved densely in real or complex arithmetic; the
    // analytic levels are 2 - 2 cos(pi (i+1) / 65) for either hopping gauge.
    for (auto const hopping : {Complex(-1.), Complex(0., 1.)}) {
        auto const result = tridiagonal(64, hopping, 0., 2.).solve(4);
        ASSERT_EQ(result.size(), 1);
        ASSERT_EQ(result[0].eigenvalues.size(), 4);
        EXPECT_EQ(result[0].method, "dense");
        for (int i = 0; i < 4; ++i)
            EXPECT_NEAR(result[0].eigenvalues[static_cast<std::size_t>(i)], 2. - 2. * std::cos(pi * (i + 1) / 65), 1e-11);
    }
    // Above the dispatch threshold both scalar types go through the Lanczos
    // driver; a linear diagonal keeps the lowest levels well separated.
    for (auto const hopping : {Complex(-.3), Complex(0., .3)}) {
        auto const p = tridiagonal(1100, hopping, .01, 2.);
        auto const full = p.solve(1100);
        auto const sparse = p.solve(4);
        ASSERT_EQ(full.size(), 1);
        ASSERT_EQ(sparse.size(), 1);
        EXPECT_EQ(full[0].method, "dense");
        EXPECT_EQ(sparse[0].method, "lanczos");
        ASSERT_EQ(sparse[0].eigenvalues.size(), 4);
        for (std::size_t i = 0; i < 4; ++i) EXPECT_NEAR(sparse[0].eigenvalues[i], full[0].eigenvalues[i], 1e-9);
    }
}

TEST(CompiledMomentum, ChargeClusteringSupportsNonHalfIntegralLocalValues) {
    auto p = ring(4, 2, .2);
    for (auto& space : p.spaces) space.quantum_numbers[0] = {.1, .3};
    auto const result = p.solve();
    for (auto const& block : result) EXPECT_NEAR(block.quantum_numbers[0], .8, 1e-15);
}

TEST(CompiledMomentum, EmptySectorReturnsNoBlocks) {
    auto p = ring(4, 5, 0.);
    EXPECT_TRUE(p.solve().empty());
}

TEST(CompiledMomentum, InputValidationRejectsSilentScientificChanges) {
    auto const original = ring(4, 1, .2);
    auto p = original;
    p.values[0] *= 2.;
    EXPECT_THROW(p.solve(), std::invalid_argument); // Not Hermitian.
    p = original;
    p.columns.push_back(0); p.values.emplace_back(1.); ++p.offsets.back();
    EXPECT_THROW(p.solve(), std::invalid_argument);
    p = original;
    p.values[0] = std::numeric_limits<double>::infinity();
    EXPECT_THROW(p.solve(), std::invalid_argument);
    p = original;
    p.columns[0] = 100;
    EXPECT_THROW(p.solve(), std::invalid_argument);
    p = original;
    p.offsets[1] = -1;
    EXPECT_THROW(p.solve(), std::invalid_argument);
    p = original;
    p.translations[0].sites[0] = p.translations[0].sites[1];
    EXPECT_THROW(p.solve(), std::invalid_argument);
    p = original;
    p.translations[0].period = 3;
    EXPECT_THROW(p.solve(), std::invalid_argument);
    p.translations[0].period = 8;
    EXPECT_THROW(p.solve(), std::invalid_argument);
    p = original;
    p.states[0] = 12;
    EXPECT_THROW(p.solve(), std::invalid_argument);
    EXPECT_THROW(original.solve(0), std::invalid_argument);
}

TEST(CompiledMomentum, TranslationBreakingPotentialIsRejected) {
    auto p = ring(4, 1, 0.);
    // Add a potential on only one product state, preserving Hermiticity.
    p.columns.insert(p.columns.begin(), 0);
    p.values.insert(p.values.begin(), 1.);
    for (std::size_t row = 1; row < p.offsets.size(); ++row) ++p.offsets[row];
    EXPECT_THROW(p.solve(), std::invalid_argument);
}

TEST(CompiledMomentum, NonconservedChargesAndIncompleteBasesAreRejected) {
    Problem p;
    p.spaces.resize(2, {{1., 1.}, {{0., 1.}}});
    p.translations.push_back({{1, 0}, 2});
    p.states = {0, 0, 0, 1, 1, 0, 1, 1};
    p.offsets = {0, 1, 1, 1, 2};
    p.columns = {3, 0};
    p.values = {1., 1.}; // |00><11| + hc changes total charge by two.
    EXPECT_THROW(p.solve(), std::invalid_argument);

    p = ring(4, 1, 0.);
    p.states.resize(12);
    p.offsets = {0, 0, 0, 0};
    p.columns.clear(); p.values.clear();
    EXPECT_THROW(p.solve(), std::invalid_argument); // Missing translated state.
    p = ring(3, 1, 0.);
    p.translations.push_back({{1, 0, 2}, 2});
    EXPECT_THROW(p.solve(), std::invalid_argument); // Noncommuting generators.
}

TEST(CompiledMomentum, DuplicateUnorderedCsrEntriesAccumulateBeforeProjection) {
    auto const original = ring(5, 2, .3);
    auto p = original;
    p.offsets = {0}; p.columns.clear(); p.values.clear();
    for (std::size_t row = 0; row + 1 < original.offsets.size(); ++row) {
        for (int copy = 0; copy < 2; ++copy)
            for (auto pos = original.offsets[row + 1]; pos-- > original.offsets[row];) {
                p.columns.push_back(original.columns[static_cast<std::size_t>(pos)]);
                p.values.push_back(.5 * original.values[static_cast<std::size_t>(pos)]);
            }
        p.offsets.push_back(static_cast<std::int64_t>(p.values.size()));
    }
    auto const expected = original.solve(), actual = p.solve();
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t block = 0; block < actual.size(); ++block) {
        ASSERT_EQ(actual[block].eigenvalues.size(), expected[block].eigenvalues.size());
        for (std::size_t i = 0; i < actual[block].eigenvalues.size(); ++i)
            EXPECT_NEAR(actual[block].eigenvalues[i], expected[block].eigenvalues[i], 1e-13);
    }
}

TEST(CompiledMomentum, ChargeOnlyResolutionPreservesTheCompleteFermionSpectrum) {
    auto p = ring(6, 2, .31);
    auto const momentum_blocks = p.solve();
    p.translations.clear();
    auto const charge_blocks = p.solve();
    ASSERT_EQ(charge_blocks.size(), 1);
    EXPECT_TRUE(charge_blocks[0].momentum.empty());
    EXPECT_EQ(charge_blocks[0].quantum_numbers, std::vector<double>{2.});
    std::vector<double> expected;
    for (auto const& block : momentum_blocks)
        expected.insert(expected.end(), block.eigenvalues.begin(), block.eigenvalues.end());
    std::sort(expected.begin(), expected.end());
    ASSERT_EQ(charge_blocks[0].eigenvalues.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) EXPECT_NEAR(charge_blocks[0].eigenvalues[i], expected[i], 1e-12);
}

TEST(CompiledMomentum, SelectedMomentaRetainTheirOriginalLabelsAndLevels) {
    auto const p = ring(6, 2, .31);
    auto const full = p.solve();
    std::vector<std::size_t> const selected{1, 4};
    auto const result = alps::ed::symmetry_spectrum(p.offsets, p.columns, p.values, p.offsets.size() - 1, p.states,
                                                    p.spaces, p.translations, 128, selected);
    ASSERT_EQ(result.size(), selected.size());
    for (std::size_t i = 0; i < selected.size(); ++i) {
        EXPECT_EQ(result[i].momentum, full[selected[i]].momentum);
        ASSERT_EQ(result[i].eigenvalues.size(), full[selected[i]].eigenvalues.size());
        for (std::size_t j = 0; j < result[i].eigenvalues.size(); ++j)
            EXPECT_NEAR(result[i].eigenvalues[j], full[selected[i]].eigenvalues[j], 1e-13);
    }
}

TEST(CompiledMomentum, SelectionControlsAndDiagnosticsCoverSparseAndDenseBlocks) {
    // A diagonal block of dimension 8 is dense for every selection; one of
    // dimension 1200 is above the Lanczos threshold for SA, LA and LM, while
    // SM (interior eigenvalues) is always solved densely.
    for (int const dimension : {8, 1200}) {
        Problem p;
        p.spaces.push_back({std::vector<double>(static_cast<std::size_t>(dimension), 1.), {}});
        for (int row = 0; row < dimension; ++row) {
            p.states.push_back(row);
            p.columns.push_back(row);
            p.values.emplace_back(row - dimension / 2 + .2);
            p.offsets.push_back(static_cast<std::int64_t>(p.values.size()));
        }
        for (auto const& which : {"SA", "LA", "SM", "LM"}) {
            SCOPED_TRACE(::testing::Message() << "dimension=" << dimension << " which=" << which);
            auto const actual = alps::ed::symmetry_spectrum(p.offsets, p.columns, p.values,
                static_cast<std::size_t>(dimension), p.states, p.spaces, p.translations, 2, {}, which, 1e-12, 2000);
            ASSERT_EQ(actual.size(), 1);
            auto const& block = actual.front();
            std::vector<double> const expected = std::string(which) == "SA"
                ? std::vector<double>{-dimension / 2 + .2, -dimension / 2 + 1.2}
                : std::string(which) == "LA" ? std::vector<double>{dimension / 2 - 1.8, dimension / 2 - .8}
                : std::string(which) == "SM" ? std::vector<double>{-.8, .2}
                : std::vector<double>{-dimension / 2 + .2, dimension / 2 - .8};
            ASSERT_EQ(block.eigenvalues.size(), 2);
            for (std::size_t i = 0; i < expected.size(); ++i) EXPECT_NEAR(block.eigenvalues[i], expected[i], 1e-10);
            EXPECT_EQ(block.nconv, 2);
            bool const lanczos = dimension == 1200 && std::string(which) != "SM";
            EXPECT_EQ(block.method, lanczos ? "lanczos" : "dense");
            if (lanczos) {
                EXPECT_DOUBLE_EQ(*block.tolerance, 1e-12);
                EXPECT_EQ(*block.max_iterations, 2000);
                EXPECT_TRUE(block.iterations.has_value()); // exact for a diagonal block: zero steps
            } else {
                EXPECT_FALSE(block.tolerance);
                EXPECT_FALSE(block.max_iterations);
                EXPECT_FALSE(block.iterations);
            }
        }
    }
}

TEST(CompiledMomentum, BlockDispatchBoundariesAndDefaultControls) {
    // One single-site block of the given dimension, through the CSR and the
    // local-operator entry points. Lanczos needs a dimension above 1000 and
    // above 8*count + 200; omitted controls are 1e-10 and min(10*n, 1000).
    struct Case {
        int dimension, count;
        bool lanczos;
        int max_iterations;
    };
    for (auto const& item : {Case{1000, 1, false, 0}, Case{1001, 1, true, 1000}, Case{1001, 101, false, 0},
                             Case{1500, 2, true, 1000}}) {
        SCOPED_TRACE(::testing::Message() << "dimension=" << item.dimension << " count=" << item.count);
        auto const dimension = item.dimension;
        Problem p;
        p.spaces.push_back({std::vector<double>(static_cast<std::size_t>(dimension), 1.), {}});
        alps::models::CompiledOperator local{{dimension}, {}};
        std::vector<Complex> matrix(static_cast<std::size_t>(dimension) * static_cast<std::size_t>(dimension));
        for (int row = 0; row < dimension; ++row) {
            p.states.push_back(row);
            for (int col : {row - 1, row, row + 1}) {
                if (col < 0 || col >= dimension) continue;
                Complex const value = col == row ? Complex(.1 * row) : Complex(.3);
                p.columns.push_back(col);
                p.values.push_back(value);
                matrix[static_cast<std::size_t>(row) * static_cast<std::size_t>(dimension) + static_cast<std::size_t>(col)] = value;
            }
            p.offsets.push_back(static_cast<std::int64_t>(p.values.size()));
        }
        local.terms.push_back({1., {{0, matrix}}, {}});
        auto const full = p.solve(dimension).front().eigenvalues;
        for (auto const& blocks : {p.solve(item.count),
                                   alps::ed::symmetry_spectrum_from_local(local, p.spaces, {}, {}, item.count)}) {
            ASSERT_EQ(blocks.size(), 1);
            auto const& block = blocks.front();
            EXPECT_EQ(block.block_dimension, static_cast<std::size_t>(dimension));
            EXPECT_EQ(block.method, item.lanczos ? "lanczos" : "dense");
            ASSERT_EQ(block.eigenvalues.size(), static_cast<std::size_t>(item.count));
            for (int i = 0; i < item.count; ++i)
                EXPECT_NEAR(block.eigenvalues[static_cast<std::size_t>(i)], full[static_cast<std::size_t>(i)], 1e-9);
            if (item.lanczos) {
                EXPECT_DOUBLE_EQ(*block.tolerance, 1e-10);
                EXPECT_EQ(*block.max_iterations, item.max_iterations);
            } else {
                EXPECT_FALSE(block.tolerance);
                EXPECT_EQ(block.nconv, item.count);
            }
        }
    }
}

TEST(CompiledMomentum, InvalidControlsAreRejectedEvenForEmptySectors) {
    auto const p = ring(4, 5, 0.);
    for (auto const tolerance : {0., -1., std::numeric_limits<double>::infinity()})
        EXPECT_THROW(alps::ed::symmetry_spectrum(p.offsets, p.columns, p.values, 0, p.states, p.spaces, p.translations,
                                                 2, {}, "SA", tolerance),
                     std::invalid_argument);
    EXPECT_THROW(alps::ed::symmetry_spectrum(p.offsets, p.columns, p.values, 0, p.states, p.spaces, p.translations, 2,
                                             {}, "XX"),
                 std::invalid_argument);
    EXPECT_THROW(alps::ed::symmetry_spectrum(p.offsets, p.columns, p.values, 0, p.states, p.spaces, p.translations, 2,
                                             {}, "SA", {}, 0),
                 std::invalid_argument);
}

TEST(CompiledMomentum, DenseSelectionRetainsMatchingEigenvectorColumns) {
    std::vector<double> const diagonal{-9., -3., -.2, .7, 2., 6.};
    for (auto const& which : {"SA", "LA", "SM", "LM"}) {
        std::vector<double> matrix(36);
        for (std::size_t i = 0; i < diagonal.size(); ++i) matrix[i * 6 + i] = diagonal[i];
        auto const result = alps::ed::dense_eigensystem(std::move(matrix), 6, true, 2, which);
        ASSERT_EQ(result.eigenvalues.size(), 2);
        ASSERT_EQ(result.eigenvectors_column_major.size(), 12);
        for (std::size_t col = 0; col < 2; ++col)
            for (std::size_t row = 0; row < 6; ++row)
                EXPECT_NEAR((diagonal[row] - result.eigenvalues[col]) * result.eigenvectors_column_major[col * 6 + row],
                            0., 1e-13);
    }
}
