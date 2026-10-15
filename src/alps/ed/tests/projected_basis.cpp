// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/tests/projected_basis.cpp.
// Adaptation: C++17 (no designated initializers, std::span or std::numbers),
// <alps/ed/...> include names; the iterative case uses a 16-site ring so the
// IETL driver's dispatch threshold selects Lanczos.
#include <alps/ed/projected_basis.hpp>
#include <alps/ed/projected_problem.hpp>
#include <alps/ed/spectrum.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <vector>

namespace {
using C = std::complex<double>;
using namespace alps::ed;
double const pi = std::acos(-1.);
EigenstateOptions options(int count, EigenstateMethod method = EigenstateMethod::automatic,
                          std::string which = "SA", std::optional<double> tolerance = {},
                          std::optional<int> max_iterations = {}) {
    EigenstateOptions o;
    o.count = count;
    o.method = method;
    o.which = std::move(which);
    o.tolerance = tolerance;
    o.max_iterations = max_iterations;
    return o;
}
// Former projected_eigenstates defaults: Lanczos requested, dense for small blocks.
EigenstateOptions lanczos(int count) { return options(count, EigenstateMethod::lanczos); }
struct Ring {
    alps::models::CompiledOperator h;
    std::vector<SpectrumLocalSpace> spaces;
    std::vector<SpectrumTranslation> translations;
    explicit Ring(int n, bool fermion = false, double phase = 0.) {
        spaces.resize(n, {{1., fermion ? -1. : 1.}, {{0., 1.}}});
        translations.resize(1);
        translations[0].period = n;
        h.site_dims.assign(n, 2);
        for (int i = 0; i < n; ++i) {
            translations[0].sites.push_back((i + 1) % n);
            int const j = (i + 1) % n;
            // c_i^dagger c_j and its adjoint in an explicit JW tensor product.
            for (int reverse = 0; reverse < 2; ++reverse) {
                int const a = reverse ? j : i, b = reverse ? i : j;
                C const coefficient = -std::polar(1., reverse ? -phase : phase) * (fermion && a > b ? -1. : 1.);
                std::vector<C> const create{0., 0., 1., 0.}, destroy{0., 1., 0., 0.};
                h.terms.push_back({coefficient,
                    {{std::min(a, b), a < b ? create : destroy}, {std::max(a, b), a < b ? destroy : create}}, {}});
                if (fermion)
                    for (int s = std::min(a, b); s < std::max(a, b); ++s) h.terms.back().jordan_wigner.push_back(s);
            }
        }
    }
    ProjectedBasis basis(double particles, std::size_t momentum) const {
        return make_projected_basis(spaces, translations, std::vector<std::optional<double>>{particles}, momentum);
    }
};
std::vector<C> dense(CooFromLocal const& coo, std::size_t rows, std::size_t cols) {
    std::vector<C> result(rows * cols);
    for (std::size_t i = 0; i < coo.data.size(); ++i)
        result[static_cast<std::size_t>(coo.rows[i]) * cols + static_cast<std::size_t>(coo.cols[i])] += coo.data[i];
    return result;
}
std::vector<C> column(std::vector<C> const& vectors, std::size_t dimension, std::size_t index) {
    return std::vector<C>(vectors.begin() + static_cast<std::ptrdiff_t>(index * dimension),
                          vectors.begin() + static_cast<std::ptrdiff_t>((index + 1) * dimension));
}
} // namespace

TEST(ProjectedBasis, ReflectionLockIncludesFermionicSignsAndCompleteSpectrum) {
    for (bool fermion : {false, true}) {
        Ring ring(6, fermion);
        std::vector<SpectrumTranslation> reflection{{{5, 4, 3, 2, 1, 0}, 2}};
        std::vector<std::optional<double>> targets{2.};
        auto parent = make_projected_basis(ring.spaces, {}, targets);
        auto oracle = prepare_projected(ring.h, parent).eigenstates(options(100)).spectrum.eigenvalues;
        std::vector<double> all;
        for (std::size_t parity = 0; parity < 2; ++parity) {
            auto basis = make_projected_basis(ring.spaces, reflection, targets, parity);
            auto problem = prepare_projected(ring.h, basis);
            auto result = problem.eigenstates(options(100));
            all.insert(all.end(), result.spectrum.eigenvalues.begin(), result.spectrum.eigenvalues.end());
            auto matrix = dense(projected_operator(ring.h, basis, basis), basis.dimension(), basis.dimension());
            for (std::size_t c = 0; c < basis.dimension(); ++c)
                for (std::size_t r = 0; r < basis.dimension(); ++r) {
                    C actual{};
                    for (std::size_t i = 0; i < basis.dimension(); ++i)
                        actual += matrix[r * basis.dimension() + i]
                            * result.eigenvectors_column_major[c * basis.dimension() + i];
                    EXPECT_NEAR(std::abs(actual - result.spectrum.eigenvalues[c]
                                             * result.eigenvectors_column_major[c * basis.dimension() + r]),
                                0., 2e-12);
                }
        }
        std::sort(all.begin(), all.end());
        ASSERT_EQ(all.size(), oracle.size());
        for (std::size_t i = 0; i < all.size(); ++i) EXPECT_NEAR(all[i], oracle[i], 2e-12);
    }
}

TEST(ProjectedBasis, OwnershipPhaseAndOneParticleDispersion) {
    Ring ring(5, true, .3);
    for (std::size_t k = 0; k < 5; ++k) {
        auto basis = ring.basis(1., k);
        ASSERT_EQ(basis.dimension(), 1);
        EXPECT_EQ(basis.parent_dimension(), 5);
        EXPECT_EQ(basis.orbit_sizes()[0], 5);
        auto const result = prepare_projected(ring.h, basis).eigenstates(lanczos(1));
        double const momentum = 2. * pi * static_cast<double>(k) / 5.;
        ASSERT_EQ(result.spectrum.eigenvalues.size(), 1);
        EXPECT_NEAR(result.spectrum.eigenvalues[0], -2. * std::cos(momentum - .3), 1e-12);
        auto const expanded = expand_projected_state(basis, std::vector<C>{1.});
        ASSERT_EQ(expanded.basis.dimension, 5);
        // Lexicographic representative is particle at site 4. Translation
        // sends it to site 0, whose coefficient is exp(-ik)/sqrt(5).
        EXPECT_NEAR(std::abs(expanded.amplitudes.front() - 1. / std::sqrt(5.)), 0., 1e-13);
        EXPECT_NEAR(std::abs(expanded.amplitudes.back() - std::polar(1. / std::sqrt(5.), -momentum)), 0., 1e-13);
    }
    auto owned = ring.basis(1., 1);
    ring.spaces[0].parity[0] = -1.;
    ring.translations[0].sites[0] = 0;
    EXPECT_EQ(owned.local_spaces()[0].parity[0], 1.);
    EXPECT_EQ(owned.translations()[0].sites[0], 1);
}

TEST(ProjectedBasis, GradedTwoParticleEnergiesAndForbiddenStabilizers) {
    Ring ring(6, true, .31);
    std::size_t total = 0;
    for (std::size_t k = 0; k < 6; ++k) {
        auto const basis = ring.basis(2., k);
        total += basis.dimension();
        auto const result = prepare_projected(ring.h, basis).eigenstates(lanczos(100));
        std::vector<double> expected;
        for (int a = 0; a < 6; ++a)
            for (int b = a + 1; b < 6; ++b)
                if ((a + b) % 6 == static_cast<int>(k))
                    expected.push_back(-2. * std::cos(2. * pi * a / 6. - .31) - 2. * std::cos(2. * pi * b / 6. - .31));
        std::sort(expected.begin(), expected.end());
        ASSERT_EQ(result.spectrum.eigenvalues.size(), expected.size());
        for (std::size_t i = 0; i < expected.size(); ++i) EXPECT_NEAR(result.spectrum.eigenvalues[i], expected[i], 1e-11);
    }
    EXPECT_EQ(total, 15);
    EXPECT_EQ(ring.basis(0., 1).dimension(), 0);
    EXPECT_EQ(ring.basis(6., 0).dimension(), 0); // cyclic permutation of 6 fermions is odd.
    EXPECT_EQ(ring.basis(6., 3).dimension(), 1);
}

TEST(ProjectedBasis, LocalCrossMomentumAndChargeChangingMeasurements) {
    Ring ring(5, true);
    auto const source = ring.basis(1., 1), target = ring.basis(1., 2);
    auto op = ring.h;
    op.terms = {{1., {{0, {0., 0., 0., 1.}}}, {}}};
    auto const matrix = dense(projected_operator(op, source, target), 1, 1);
    C const expected = std::polar(.2, 2. * pi / 5.);
    EXPECT_NEAR(std::abs(matrix[0] - expected), 0., 1e-13);
    EXPECT_NEAR(std::abs(apply_projected_operator(op, source, target, std::vector<C>{1.})[0] - expected), 0., 1e-13);
    EXPECT_NEAR(std::abs(expect_projected_operator(op, source, std::vector<C>{1.}) - .2), 0., 1e-13);
    op.terms[0].factors[0].matrix = {0., 1., 0., 0.};
    auto const vacuum = ring.basis(0., 0);
    auto const removed = apply_projected_operator(op, source, vacuum, std::vector<C>{1.});
    ASSERT_EQ(removed.size(), 1);
    EXPECT_NEAR(std::abs(removed[0] - std::polar(1. / std::sqrt(5.), -2. * pi / 5.)), 0., 1e-13);
}

TEST(ProjectedBasis, GeneralCancellationAndInvalidHamiltonians) {
    Ring ring(4);
    auto const basis = ring.basis(2., 0);
    auto const reference = prepare_projected(ring.h, basis).eigenstates(lanczos(100)).spectrum.eigenvalues;
    // A structurally different cancellation must fall back, not be rejected.
    ring.h.terms.insert(ring.h.terms.end(), {{1., {{0, {0., 0., 0., 1.}}}, {}}, {-.5, {{0, {0., 0., 0., 2.}}}, {}}});
    auto const actual = prepare_projected(ring.h, basis).eigenstates(lanczos(100)).spectrum.eigenvalues;
    ASSERT_EQ(actual.size(), reference.size());
    for (std::size_t i = 0; i < actual.size(); ++i) EXPECT_NEAR(actual[i], reference[i], 1e-12);
    ring.h.terms.back().coefficient = -.4;
    EXPECT_THROW(prepare_projected(ring.h, basis).eigenstates(lanczos(1)), std::invalid_argument);
    ring.h.terms.back().coefficient = std::numeric_limits<double>::infinity();
    EXPECT_THROW(prepare_projected(ring.h, basis).eigenstates(lanczos(1)), std::invalid_argument);
    Ring broken(4);
    broken.h.terms[0].coefficient *= C(0., 1.);
    EXPECT_THROW(prepare_projected(broken.h, basis).eigenstates(lanczos(1)), std::invalid_argument);
}

TEST(ProjectedBasis, HugeProductSpaceAndEmptySector) {
    Ring ring(70, true);
    auto const basis = ring.basis(0., 0);
    EXPECT_EQ(basis.parent_dimension(), 1);
    EXPECT_EQ(basis.dimension(), 1);
    EXPECT_EQ(prepare_projected(ring.h, basis).eigenstates(lanczos(1)).spectrum.eigenvalues, std::vector<double>{0.});
    auto const empty = ring.basis(71., 0);
    EXPECT_EQ(empty.dimension(), 0);
    EXPECT_TRUE(prepare_projected(ring.h, empty).eigenstates(lanczos(1)).spectrum.eigenvalues.empty());
    EXPECT_THROW(prepare_projected(ring.h, empty).eigenstates(lanczos(0)), std::invalid_argument);
    auto invalid_tolerance = lanczos(1);
    invalid_tolerance.tolerance = 0.;
    EXPECT_THROW(prepare_projected(ring.h, empty).eigenstates(invalid_tolerance), std::invalid_argument);
}

TEST(ProjectedBasis, JointSpectrumMatchesIndependentParentProjection) {
    for (bool fermion : {false, true}) {
        Ring ring(6, fermion, .27);
        auto const parent = enumerate_product_basis(ring.h.site_dims);
        std::vector<std::int64_t> const states(parent.states.begin(), parent.states.end());
        std::vector<std::vector<double>> parity;
        for (auto const& space : ring.spaces) parity.push_back(space.parity);
        auto const coo = build_coo_from_local(ring.h, states, parent.dimension, parity);
        std::vector<std::map<std::int64_t, C>> rows(parent.dimension);
        for (std::size_t i = 0; i < coo.data.size(); ++i)
            rows[static_cast<std::size_t>(coo.rows[i])][coo.cols[i]] += coo.data[i];
        std::vector<std::int64_t> offsets{0}, columns;
        std::vector<C> values;
        for (auto const& row : rows) {
            for (auto const& entry : row) {
                columns.push_back(entry.first);
                values.push_back(entry.second);
            }
            offsets.push_back(static_cast<std::int64_t>(values.size()));
        }
        auto const expected = symmetry_spectrum(offsets, columns, values, parent.dimension, states, ring.spaces,
                                                ring.translations, 100);
        auto const actual = symmetry_spectrum_from_local(ring.h, ring.spaces, ring.translations, {}, 100);
        ASSERT_EQ(actual.size(), expected.size());
        for (std::size_t b = 0; b < actual.size(); ++b) {
            EXPECT_EQ(actual[b].momentum, expected[b].momentum);
            EXPECT_EQ(actual[b].quantum_numbers, expected[b].quantum_numbers);
            ASSERT_EQ(actual[b].eigenvalues.size(), expected[b].eigenvalues.size());
            for (std::size_t i = 0; i < actual[b].eigenvalues.size(); ++i)
                EXPECT_NEAR(actual[b].eigenvalues[i], expected[b].eigenvalues[i], 1e-11);
        }
    }
}

TEST(ProjectedBasis, IterativeVectorsAndSelectionRetainTheirResiduals) {
    // Half filling of a 16-site ring at momentum zero: the block (about 810
    // states) exceeds the Lanczos dispatch threshold of the IETL driver.
    Ring ring(16);
    auto const basis = ring.basis(8., 0);
    ASSERT_GT(basis.dimension(), 8 * 3 + 200);
    for (auto const* which : {"SA", "LA"}) {
        auto const problem = prepare_projected(ring.h, basis);
        auto const sparse = problem.eigenstates(options(3, EigenstateMethod::lanczos, which, 1e-12, 3000));
        auto const dense_result = problem.eigenstates(options(3, EigenstateMethod::dense, which));
        ASSERT_EQ(sparse.spectrum.method, "lanczos");
        ASSERT_EQ(sparse.residual_norms.size(), 3);
        for (std::size_t i = 0; i < 3; ++i) {
            EXPECT_NEAR(sparse.spectrum.eigenvalues[i], dense_result.spectrum.eigenvalues[i], 1e-10);
            EXPECT_LT(sparse.residual_norms[i], 1e-9);
            auto const vector = column(sparse.eigenvectors_column_major, basis.dimension(), i);
            EXPECT_NEAR(expect_projected_operator(ring.h, basis, vector).real(), sparse.spectrum.eigenvalues[i], 1e-10);
        }
    }
}

TEST(ProjectedBasis, TwoDirectionsHavePositiveTranslationEigenphases) {
    std::vector<SpectrumLocalSpace> const spaces(6, {{1., 1.}, {{0., 1.}}});
    std::vector<SpectrumTranslation> const translations{{{2, 3, 4, 5, 0, 1}, 3}, {{1, 0, 3, 2, 5, 4}, 2}};
    for (std::size_t k = 0; k < 6; ++k) {
        auto const basis = make_projected_basis(spaces, translations, std::vector<std::optional<double>>{1.}, k);
        ASSERT_EQ(basis.dimension(), 1);
        auto const expanded = expand_projected_state(basis, std::vector<C>{1.});
        for (std::size_t direction = 0; direction < 2; ++direction) {
            std::vector<C> shifted(6);
            // One-particle lexicographic rows enumerate sites in reverse.
            for (std::size_t source = 0; source < 6; ++source)
                shifted[5 - translations[direction].sites[source]] = expanded.amplitudes[5 - source];
            for (std::size_t r = 0; r < 6; ++r)
                EXPECT_NEAR(std::abs(shifted[r] - std::polar(1., basis.momentum()[direction]) * expanded.amplitudes[r]),
                            0., 1e-13);
        }
    }
}

TEST(ProjectedBasis, ValidationIncludesUnfilteredChargeChangingTermsAndDimensions) {
    Ring ring(4);
    auto const basis = ring.basis(2., 0);
    auto broken = ring.h;
    broken.terms.push_back({1., {{0, {0., 1., 1., 0.}}}, {}});
    EXPECT_THROW(prepare_projected(broken, basis).eigenstates(lanczos(1)), std::invalid_argument);
    EXPECT_THROW(make_projected_basis(ring.spaces, ring.translations, {}, 4), std::invalid_argument);
    EXPECT_THROW(make_projected_basis(ring.spaces, ring.translations, std::vector<std::optional<double>>{1., 2.}),
                 std::invalid_argument);
    auto badtranslations = ring.translations;
    badtranslations[0].sites[0] = badtranslations[0].sites[1];
    EXPECT_THROW(make_projected_basis(ring.spaces, badtranslations), std::invalid_argument);
    EXPECT_THROW(expect_projected_operator(ring.h, basis, std::vector<C>{}), std::invalid_argument);
    auto other = ring.spaces;
    other[0].parity[0] = -1.;
    auto const incompatible = make_projected_basis(other, {}, std::vector<std::optional<double>>{2.});
    EXPECT_THROW(projected_operator(ring.h, basis, incompatible), std::invalid_argument);
}

TEST(ProjectedBasis, FixedChargeAnchorsRejectNearToleranceCrossBlockCoupling) {
    std::vector<SpectrumLocalSpace> const spaces{{{1., 1., 1.}, {{0., .75e-9, 1.5e-9}}}};
    alps::models::CompiledOperator h;
    h.site_dims = {3};
    h.terms = {{1., {{0, {0., 0., 0., 0., 0., 1., 0., 1., 0.}}}, {}}};
    auto const basis = make_projected_basis(spaces, {});
    EXPECT_THROW(prepare_projected(h, basis).eigenstates(lanczos(3)), std::invalid_argument);
    EXPECT_THROW(symmetry_spectrum_from_local(h, spaces, {}, {}, 3), std::invalid_argument);
}
