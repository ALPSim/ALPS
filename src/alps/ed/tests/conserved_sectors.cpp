// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/tests/conserved_sectors.cpp.
// Adaptation: C++17 (no designated initializers), <alps/ed/...> include names.
#include <alps/ed/projected_basis.hpp>
#include <alps/ed/projected_problem.hpp>
#include <alps/ed/spectrum.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
using C = std::complex<double>;
using Input = alps::models::CompiledOperator;
using namespace alps::ed;
// Former projected_eigenstates defaults: Lanczos requested, dense for small blocks.
EigenstateOptions lanczos(int count) {
    EigenstateOptions o;
    o.count = count;
    o.method = EigenstateMethod::lanczos;
    return o;
}

std::vector<C> const identity{1., 0., 0., 1.};
std::vector<C> const sx{0., .5, .5, 0.};
std::vector<C> const sy{0., C(0., -.5), C(0., .5), 0.};
std::vector<C> const sz{.5, 0., 0., -.5};
std::vector<C> const pup{1., 0., 0., 0.};
std::vector<C> const pdown{0., 0., 0., 1.};
std::vector<SpectrumLocalSpace> const spaces(2, {{1., 1.}, {{.5, -.5}}});
std::vector<std::optional<double>> const zero_charge{0.};

std::vector<SpectrumTranslation> translations(bool momentum) {
    return momentum ? std::vector<SpectrumTranslation>{{{1, 0}, 2}} : std::vector<SpectrumTranslation>{};
}

Input empty_input() {
    Input input;
    input.site_dims = {2, 2};
    return input;
}

void term(Input& input, C coefficient, std::vector<int> sites, std::vector<std::vector<C>> operators) {
    alps::models::LocalTerm result{coefficient, {}, {}};
    for (std::size_t i = 0; i < sites.size(); ++i) result.factors.push_back({sites[i], std::move(operators.at(i))});
    input.terms.push_back(std::move(result));
}

Input transverse(double offset, double factor_scale = 1.) {
    auto input = empty_input();
    term(input, offset, {0}, {identity});
    auto scaled = sx;
    for (auto& entry : scaled) entry *= factor_scale;
    term(input, 1. / factor_scale, {0}, {scaled});
    term(input, 1. / factor_scale, {1}, {scaled});
    return input;
}

// Encode independently specified 4x4 matrix entries; do not use the local
// operator compiler to construct the CSR reference. Product rows are
// |up up>, |up down>, |down up>, |down down>.
struct MatrixInput {
    std::vector<std::int64_t> offsets{0}, columns;
    std::vector<C> values;
    std::vector<std::int64_t> const states{0, 0, 0, 1, 1, 0, 1, 1};

    explicit MatrixInput(std::array<C, 16> const& matrix) {
        for (int row = 0; row < 4; ++row) {
            for (int col = 0; col < 4; ++col)
                if (matrix[static_cast<std::size_t>(4 * row + col)] != C{}) {
                    columns.push_back(col);
                    values.push_back(matrix[static_cast<std::size_t>(4 * row + col)]);
                }
            offsets.push_back(static_cast<std::int64_t>(values.size()));
        }
    }

    std::vector<SpectrumBlock> solve(bool momentum) const {
        return symmetry_spectrum(offsets, columns, values, 4, states, spaces, translations(momentum), 4);
    }
};

Input transverse_contributions(std::vector<double> const& coefficients, bool imaginary) {
    auto input = empty_input();
    auto const& flip = imaginary ? sy : sx;
    for (double coefficient : coefficients) {
        term(input, coefficient, {0}, {flip});
        term(input, coefficient, {1}, {flip});
    }
    return input;
}

MatrixInput duplicate_transverse(std::vector<double> const& coefficients, bool imaginary, bool diagonal = false,
                                 bool reverse_odd_rows = false) {
    MatrixInput input(std::array<C, 16>{});
    input.offsets = {0};
    for (int row = 0; row < 4; ++row) {
        if (diagonal) {
            input.columns.push_back(row);
            input.values.push_back(row == 0 || row == 3 ? .25 : -.25);
        }
        for (int col = 0; col < 4; ++col) {
            if ((row ^ col) != 1 && (row ^ col) != 2) continue;
            // Each pair differs by one spin. Sy has -i/2 above the
            // diagonal and +i/2 below; Sx has +1/2 in both positions.
            C const entry = imaginary ? C(0., row < col ? -.5 : .5) : C(.5);
            for (std::size_t i = 0; i < coefficients.size(); ++i) {
                auto const index = reverse_odd_rows && row % 2 ? coefficients.size() - 1 - i : i;
                input.columns.push_back(col);
                input.values.push_back(coefficients[index] * entry);
            }
        }
        input.offsets.push_back(static_cast<std::int64_t>(input.values.size()));
    }
    return input;
}

std::vector<double> energies(std::vector<SpectrumBlock> const& blocks) {
    std::vector<double> result;
    for (auto const& block : blocks) result.insert(result.end(), block.eigenvalues.begin(), block.eigenvalues.end());
    std::sort(result.begin(), result.end());
    return result;
}

void expect_energies(std::vector<double> const& actual, std::vector<double> const& expected) {
    ASSERT_EQ(actual.size(), expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i) EXPECT_NEAR(actual[i], expected[i], 1e-13);
}
} // namespace

TEST(ConservedSectors, LocalChargeEscapeCannotBeHiddenByIdentityOrFactorScaling) {
    for (bool momentum : {false, true})
        for (double offset : {0., 0x1p50})
            for (double factor_scale : {1., 0x1p-40}) {
                SCOPED_TRACE(::testing::Message() << "momentum=" << momentum << " offset=" << offset
                                                  << " factor_scale=" << factor_scale);
                auto const input = transverse(offset, factor_scale);
                auto const lattice = translations(momentum);
                auto const basis = make_projected_basis(spaces, lattice, zero_charge);
                // [Q, offset*I+Sx0+Sx1] = i(Sy0+Sy1). Every escape from
                // Sz=0 has amplitude 1/2, even with a huge energy offset.
                EXPECT_THROW(prepare_projected(input, basis).eigenstates(lanczos(2)), std::invalid_argument);
                EXPECT_THROW(symmetry_spectrum_from_local(input, spaces, lattice, zero_charge, 2), std::invalid_argument);
            }
}

TEST(ConservedSectors, CsrChargeEscapeCannotBeHiddenByIdentity) {
    for (bool momentum : {false, true})
        for (double offset : {0., 0x1p50}) {
            SCOPED_TRACE(::testing::Message() << "momentum=" << momentum << " offset=" << offset);
            std::array<C, 16> matrix{};
            for (int row = 0; row < 4; ++row)
                for (int col = 0; col < 4; ++col)
                    matrix[static_cast<std::size_t>(4 * row + col)]
                        = row == col ? offset : ((row ^ col) == 1 || (row ^ col) == 2) ? .5 : 0.;
            // The explicit CSR has the same eight forbidden 1/2 entries;
            // no local decomposition or selected-basis assembly is used.
            EXPECT_THROW(MatrixInput(matrix).solve(momentum), std::invalid_argument);
        }
}

TEST(ConservedSectors, LowLevelProjectedSolveRetainsAnExplicitInvariantParentSector) {
    auto input = empty_input();
    term(input, 1., {0, 1}, {sx, sx});
    term(input, -1., {0, 1}, {sy, sy});
    // (S+S+ + S-S-)/2 breaks global Sz but annihilates the selected Sz=0
    // parent sector. This lower-level API promises selected-sector
    // invariance, unlike a model declaration of full-operator conservation.
    for (bool momentum : {false, true}) {
        auto const basis = make_projected_basis(spaces, translations(momentum), zero_charge);
        auto const result = prepare_projected(input, basis).eigenstates(lanczos(2));
        expect_energies(result.spectrum.eigenvalues, momentum ? std::vector<double>{0.} : std::vector<double>{0., 0.});
    }
}

TEST(ConservedSectors, CartesianCancellationRetainsAnalyticChargeAndMomentumSpectra) {
    auto input = empty_input();
    term(input, 1., {0, 1}, {sx, sx});
    term(input, 1., {0, 1}, {sy, sy});
    // SxSx+SySy cancels both pair-creation entries. In Sz=0 the only
    // entries are H[up down,down up]=H[down up,up down]=1/2.
    std::array<C, 16> matrix{};
    matrix[1 * 4 + 2] = matrix[2 * 4 + 1] = .5;
    for (bool momentum : {false, true}) {
        auto const lattice = translations(momentum);
        auto const basis = make_projected_basis(spaces, lattice, zero_charge);
        expect_energies(prepare_projected(input, basis).eigenstates(lanczos(2)).spectrum.eigenvalues,
                        momentum ? std::vector<double>{.5} : std::vector<double>{-.5, .5});
        expect_energies(energies(symmetry_spectrum_from_local(input, spaces, lattice, zero_charge, 2)), {-.5, .5});
        expect_energies(energies(MatrixInput(matrix).solve(momentum)), {-.5, 0., 0., .5});
    }
}

TEST(ConservedSectors, CrossSupportCancellationIsNotRejectedTermByTerm) {
    for (bool reverse : {false, true})
        for (bool momentum : {false, true}) {
            SCOPED_TRACE(::testing::Message() << "reverse=" << reverse << " momentum=" << momentum);
            auto input = empty_input();
            term(input, 1., {0, 1}, {sz, sz});
            auto cancellation = empty_input();
            term(cancellation, 0x1p50, {0, 1}, {sx, pup});
            term(cancellation, 0x1p50, {0, 1}, {sx, pdown});
            term(cancellation, -0x1p50, {0}, {sx});
            for (int i = 0; i < 3; ++i) {
                int const index = reverse ? 2 - i : i;
                input.terms.push_back(cancellation.terms[static_cast<std::size_t>(index)]);
            }
            // Pup+Pdown=I removes all large, individually charge-changing
            // terms. Remaining Sz0*Sz1 has energy -1/4 throughout Sz=0.
            auto const lattice = translations(momentum);
            auto const basis = make_projected_basis(spaces, lattice, zero_charge);
            expect_energies(prepare_projected(input, basis).eigenstates(lanczos(2)).spectrum.eigenvalues,
                            momentum ? std::vector<double>{-.25} : std::vector<double>{-.25, -.25});
            expect_energies(energies(symmetry_spectrum_from_local(input, spaces, lattice, zero_charge, 2)), {-.25, -.25});
        }
}

TEST(ConservedSectors, ExplicitChargeChangingProjectionRemainsAvailable) {
    auto const input = transverse(0.);
    std::vector<std::optional<double>> const plus_charge{1.};
    for (bool momentum : {false, true}) {
        auto const lattice = translations(momentum);
        auto const source = make_projected_basis(spaces, lattice, zero_charge);
        auto const target = make_projected_basis(spaces, lattice, plus_charge);
        std::vector<C> const vector(source.dimension(), 1. / std::sqrt(static_cast<double>(source.dimension())));
        auto const projected = apply_projected_operator(input, source, source, vector);
        ASSERT_EQ(projected.size(), source.dimension());
        for (auto amplitude : projected) EXPECT_EQ(amplitude, C{});
        // [1/2,1/2] acting on the normalized symmetric product vector,
        // or [1/sqrt(2)] in momentum zero, yields the same amplitude.
        auto const changed = apply_projected_operator(input, source, target, vector);
        ASSERT_EQ(changed.size(), 1);
        EXPECT_NEAR(std::abs(changed[0] - 1. / std::sqrt(2.)), 0., 1e-13);
        EXPECT_EQ(expect_projected_operator(input, source, vector), C{});
    }
}

TEST(ConservedSectors, LocalEscapeValidationRetainsAResidualLostByDoubleSummation) {
    for (bool momentum : {false, true})
        for (bool imaginary : {false, true})
            for (double scale : {1., 0x1p-500})
                for (bool cancel_large_first : {false, true}) {
                    SCOPED_TRACE(::testing::Message() << "momentum=" << momentum << " imaginary=" << imaginary
                                                      << " scale=" << scale << " cancel_large_first=" << cancel_large_first);
                    auto const coefficients = cancel_large_first
                        ? std::vector<double>{0x1p54 * scale, -0x1p54 * scale, scale}
                        : std::vector<double>{0x1p54 * scale, scale, -0x1p54 * scale};
                    auto const input = transverse_contributions(coefficients, imaginary);
                    auto const lattice = translations(momentum);
                    auto const basis = make_projected_basis(spaces, lattice, zero_charge);
                    // Binary-rational coefficient sum is exactly scale.
                    // A forbidden amplitude is scale/2 (real or imaginary),
                    // although left-to-right double addition can produce zero.
                    EXPECT_THROW(prepare_projected(input, basis).eigenstates(lanczos(2)), std::invalid_argument);
                    EXPECT_THROW(symmetry_spectrum_from_local(input, spaces, lattice, zero_charge, 2),
                                 std::invalid_argument);
                }
}

TEST(ConservedSectors, CsrDuplicateEscapeValidationRetainsAResidualLostByDoubleSummation) {
    for (bool momentum : {false, true})
        for (bool imaginary : {false, true})
            for (double scale : {1., 0x1p-500})
                for (bool cancel_large_first : {false, true}) {
                    SCOPED_TRACE(::testing::Message() << "momentum=" << momentum << " imaginary=" << imaginary
                                                      << " scale=" << scale << " cancel_large_first=" << cancel_large_first);
                    auto const coefficients = cancel_large_first
                        ? std::vector<double>{0x1p54 * scale, -0x1p54 * scale, scale}
                        : std::vector<double>{0x1p54 * scale, scale, -0x1p54 * scale};
                    // The CSR contains three explicit values at each
                    // forbidden row/column pair. Their exact sum is scale/2;
                    // validation must precede a lossy duplicate collapse.
                    EXPECT_THROW(duplicate_transverse(coefficients, imaginary).solve(momentum), std::invalid_argument);
                }
}

TEST(ConservedSectors, LocalExactlyCancellingEscapesRemainValidInEveryTermOrder) {
    for (bool momentum : {false, true})
        for (bool imaginary : {false, true})
            for (bool cancel_large_first : {false, true}) {
                SCOPED_TRACE(::testing::Message() << "momentum=" << momentum << " imaginary=" << imaginary
                                                  << " cancel_large_first=" << cancel_large_first);
                auto const coefficients = cancel_large_first ? std::vector<double>{0x1p54, -0x1p54, 1., -1.}
                                                             : std::vector<double>{0x1p54, 1., -0x1p54, -1.};
                auto input = transverse_contributions(coefficients, imaginary);
                term(input, 1., {0, 1}, {sz, sz});
                auto const lattice = translations(momentum);
                auto const basis = make_projected_basis(spaces, lattice, zero_charge);
                // A+1-A-1 is exactly zero. Rejecting the rounded residual
                // would reject a conserving Hamiltonian; only Sz0*Sz1 remains.
                auto const expected = momentum ? std::vector<double>{-.25} : std::vector<double>{-.25, -.25};
                EXPECT_NO_THROW(expect_energies(prepare_projected(input, basis).eigenstates(lanczos(2)).spectrum.eigenvalues,
                                                expected));
                EXPECT_NO_THROW(expect_energies(energies(symmetry_spectrum_from_local(input, spaces, lattice, zero_charge, 2)),
                                                (std::vector<double>{-.25, -.25})));
            }
}

TEST(ConservedSectors, LocalExactCancellationPreservesAnUnrestrictedChargeParent) {
    for (bool momentum : {false, true})
        for (bool imaginary : {false, true})
            for (bool cancel_large_first : {false, true}) {
                SCOPED_TRACE(::testing::Message() << "momentum=" << momentum << " imaginary=" << imaginary
                                                  << " cancel_large_first=" << cancel_large_first);
                auto const coefficients = cancel_large_first ? std::vector<double>{0x1p54, -0x1p54, 1., -1.}
                                                             : std::vector<double>{0x1p54, 1., -0x1p54, -1.};
                auto input = transverse_contributions(coefficients, imaginary);
                term(input, 1., {0, 1}, {sz, sz});
                auto const lattice = translations(momentum);
                auto const basis = make_projected_basis(spaces, lattice, {});
                // Keep all declared charge sectors in the parent: a rounded
                // forbidden residual cannot disappear by selecting Sz=0.
                // The exact operator is Sz0*Sz1, diagonal with entries
                // (1/4,-1/4,-1/4,1/4). Momentum zero keeps |up up>,
                // (|up down>+|down up>)/sqrt(2), and |down down>.
                auto const expected = momentum ? std::vector<double>{-.25, .25, .25}
                                              : std::vector<double>{-.25, -.25, .25, .25};
                ASSERT_EQ(basis.dimension(), expected.size());
                EXPECT_NO_THROW(expect_energies(prepare_projected(input, basis).eigenstates(lanczos(4)).spectrum.eigenvalues,
                                                expected));
                // Charge-resolved spectra retain the same four physical
                // levels when collecting all requested momentum sectors.
                EXPECT_NO_THROW(expect_energies(energies(symmetry_spectrum_from_local(input, spaces, lattice, {}, 4)),
                                                (std::vector<double>{-.25, -.25, .25, .25})));
            }
}

TEST(ConservedSectors, CsrExactlyCancellingDuplicatesRemainValidInEveryEntryOrder) {
    for (bool momentum : {false, true})
        for (bool imaginary : {false, true})
            for (bool cancel_large_first : {false, true})
                for (bool reverse_odd_rows : {false, true}) {
                    SCOPED_TRACE(::testing::Message() << "momentum=" << momentum << " imaginary=" << imaginary
                                                      << " cancel_large_first=" << cancel_large_first
                                                      << " reverse_odd_rows=" << reverse_odd_rows);
                    auto const coefficients = cancel_large_first ? std::vector<double>{0x1p54, -0x1p54, 1., -1.}
                                                                 : std::vector<double>{0x1p54, 1., -0x1p54, -1.};
                    // All forbidden duplicate entries cancel exactly. Each
                    // row may order those entries independently; differing
                    // rounded residues must not manufacture non-Hermiticity.
                    // The remaining diagonal is Sz0*Sz1 on four states.
                    EXPECT_NO_THROW(expect_energies(
                        energies(duplicate_transverse(coefficients, imaginary, true, reverse_odd_rows).solve(momentum)),
                        (std::vector<double>{-.25, -.25, .25, .25})));
                }
}

TEST(ConservedSectors, LocalTinyEscapesAreNotAnAmplitudeTolerance) {
    for (bool momentum : {false, true})
        for (bool imaginary : {false, true})
            for (bool scale_in_factor : {false, true}) {
                SCOPED_TRACE(::testing::Message() << "momentum=" << momentum << " imaginary=" << imaginary
                                                  << " scale_in_factor=" << scale_in_factor);
                auto input = empty_input();
                term(input, 1., {0}, {identity});
                auto flip = imaginary ? sy : sx;
                if (scale_in_factor)
                    for (auto& entry : flip) entry *= 0x1p-1000;
                double const coefficient = scale_in_factor ? 0x1p500 : 0x1p-500;
                term(input, coefficient, {0}, {flip});
                term(input, coefficient, {1}, {flip});
                auto const lattice = translations(momentum);
                auto const basis = make_projected_basis(spaces, lattice, zero_charge);
                // All supplied numbers and the final amplitude 2^-501 are
                // representable binary doubles. Moving the scale into a
                // local factor or adding I cannot change charge conservation.
                EXPECT_THROW(prepare_projected(input, basis).eigenstates(lanczos(2)), std::invalid_argument);
                EXPECT_THROW(symmetry_spectrum_from_local(input, spaces, lattice, zero_charge, 2), std::invalid_argument);
            }
}

TEST(ConservedSectors, CsrTinyEscapesAreNotAnAmplitudeTolerance) {
    for (bool momentum : {false, true})
        for (bool imaginary : {false, true}) {
            SCOPED_TRACE(::testing::Message() << "momentum=" << momentum << " imaginary=" << imaginary);
            // Explicit nonzero cross-charge entries 2^-501 accompany a
            // finite, translation-invariant Sz0*Sz1 diagonal. No duplicates
            // or floating-point cancellation are needed in this case.
            EXPECT_THROW(duplicate_transverse({0x1p-500}, imaginary, true).solve(momentum), std::invalid_argument);
        }
}
