// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Lock for alps::ed::lanczos_csr, the deflated IETL Lanczos driver. Checks
// after skilledwolf/alpsx src/solvers/diag/tests/{lanczos_csr_lock,
// lanczos_complex_csr}.cpp (Spectra-based), restated for this driver's
// contract: eigenvalue correctness against analytic tridiagonals, the
// column-major eigenvector layout via residuals, orthonormality, degenerate
// multiplets through deflation, selection names and argument validation.

#include <alps/ed/eigensolver.hpp>
#include <alps/ed/lanczos.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using C = std::complex<double>;
double const pi = std::acos(-1.);
double conjugate(double z) { return z; }
C conjugate(C z) { return std::conj(z); }

template<class Scalar> struct Csr {
    std::vector<std::int64_t> offsets{0}, indices;
    std::vector<Scalar> values;
    std::int64_t dim = 0;
    void add(std::int64_t col, Scalar value) { indices.push_back(col); values.push_back(value); }
    void end_row() { offsets.push_back(static_cast<std::int64_t>(values.size())); ++dim; }
};

// Open chain: diagonal `shift`, hopping t between neighbours; eigenvalues
// shift - 2|t| cos(j pi/(n+1)), j = 1..n.
template<class Scalar> Csr<Scalar> chain(int n, Scalar t, double shift = 0.) {
    Csr<Scalar> m;
    for (int row = 0; row < n; ++row) {
        if (row > 0) m.add(row - 1, t);
        if (shift != 0.) m.add(row, Scalar(shift));
        if (row + 1 < n) m.add(row + 1, conjugate(t));
        m.end_row();
    }
    return m;
}

template<class Scalar> std::vector<double> chain_spectrum(int n, Scalar t, double shift = 0.) {
    std::vector<double> out;
    for (int j = 1; j <= n; ++j) out.push_back(shift - 2. * std::abs(t) * std::cos(j * pi / (n + 1)));
    std::sort(out.begin(), out.end());
    return out;
}

template<class Scalar>
void check_eigenpairs(Csr<Scalar> const& m, alps::ed::LanczosResult<Scalar> const& r,
                      double residual_tolerance = 1e-7) {
    auto const n = static_cast<std::size_t>(m.dim), k = r.eigenvalues.size();
    ASSERT_EQ(r.nconv, static_cast<int>(k));
    ASSERT_TRUE(std::is_sorted(r.eigenvalues.begin(), r.eigenvalues.end()));
    ASSERT_EQ(r.eigenvectors.size(), n * k);
    ASSERT_EQ(r.residual_norms.size(), k);
    for (std::size_t col = 0; col < k; ++col) {
        double residual = 0.;
        for (std::size_t row = 0; row < n; ++row) {
            Scalar mv{};
            for (auto p = m.offsets[row]; p < m.offsets[row + 1]; ++p)
                mv += m.values[static_cast<std::size_t>(p)]
                    * r.eigenvectors[col * n + static_cast<std::size_t>(m.indices[static_cast<std::size_t>(p)])];
            residual += std::norm(mv - r.eigenvalues[col] * r.eigenvectors[col * n + row]);
        }
        residual = std::sqrt(residual);
        EXPECT_LT(residual, residual_tolerance) << "column " << col;
        EXPECT_NEAR(r.residual_norms[col], residual, 1e-12) << "column " << col;
        for (std::size_t other = 0; other < k; ++other) {
            Scalar overlap{};
            for (std::size_t row = 0; row < n; ++row)
                overlap += conjugate(r.eigenvectors[col * n + row]) * r.eigenvectors[other * n + row];
            EXPECT_NEAR(std::abs(overlap - Scalar(col == other ? 1. : 0.)), 0., 1e-9)
                << "columns " << col << " and " << other;
        }
    }
}
} // namespace

TEST(LanczosCsr, LowestEigenpairsOfRealAndComplexChains) {
    {
        auto const m = chain<double>(200, -1.);
        auto const r = alps::ed::lanczos_csr(m.offsets, m.indices, m.values, m.dim, 3);
        EXPECT_EQ(r.tolerance, 1e-10);
        EXPECT_EQ(r.max_iterations, 1000);
        EXPECT_GT(r.iterations, 0);
        EXPECT_GE(r.runs, 2);   // one harvest and at least one verification
        auto const expected = chain_spectrum(200, -1.);
        ASSERT_EQ(r.eigenvalues.size(), 3);
        for (std::size_t i = 0; i < 3; ++i) EXPECT_NEAR(r.eigenvalues[i], expected[i], 1e-8);
        check_eigenpairs(m, r);
        auto const values_only = alps::ed::lanczos_csr(m.offsets, m.indices, m.values, m.dim, 3, "SA", 0., 0, false);
        EXPECT_EQ(values_only.eigenvalues, r.eigenvalues);
        EXPECT_TRUE(values_only.eigenvectors.empty());
        EXPECT_TRUE(values_only.residual_norms.empty());
    }
    {
        C const t = -std::polar(1., 0.41);
        auto const m = chain<C>(150, t, 0.125);
        auto const r = alps::ed::lanczos_csr(m.offsets, m.indices, m.values, m.dim, 4, "SA", 1e-12, 2000);
        EXPECT_EQ(r.tolerance, 1e-12);
        EXPECT_EQ(r.max_iterations, 2000);
        auto const expected = chain_spectrum(150, t, 0.125);
        ASSERT_EQ(r.eigenvalues.size(), 4);
        for (std::size_t i = 0; i < 4; ++i) EXPECT_NEAR(r.eigenvalues[i], expected[i], 1e-8);
        check_eigenpairs(m, r);
    }
}

TEST(LanczosCsr, LargestAndLargestMagnitudeSelection) {
    auto const m = chain<double>(200, -1., 0.5);   // spectrum in (-1.5, 2.5)
    auto const expected = chain_spectrum(200, -1., 0.5);
    auto const top = alps::ed::lanczos_csr(m.offsets, m.indices, m.values, m.dim, 3, "LA");
    ASSERT_EQ(top.eigenvalues.size(), 3);
    for (std::size_t i = 0; i < 3; ++i) EXPECT_NEAR(top.eigenvalues[i], expected[197 + i], 1e-8);
    check_eigenpairs(m, top);
    auto const magnitude = alps::ed::lanczos_csr(m.offsets, m.indices, m.values, m.dim, 3, "LM");
    ASSERT_EQ(magnitude.eigenvalues.size(), 3);
    for (std::size_t i = 0; i < 3; ++i) EXPECT_NEAR(magnitude.eigenvalues[i], expected[197 + i], 1e-8);
    check_eigenpairs(m, magnitude);
    // A symmetric spectrum: LM takes from both ends.
    auto const symmetric = chain<double>(200, -1.);
    auto const both = alps::ed::lanczos_csr(symmetric.offsets, symmetric.indices, symmetric.values, symmetric.dim, 2, "LM");
    ASSERT_EQ(both.eigenvalues.size(), 2);
    auto const reference = chain_spectrum(200, -1.);
    EXPECT_NEAR(both.eigenvalues[0], reference.front(), 1e-8);
    EXPECT_NEAR(both.eigenvalues[1], reference.back(), 1e-8);
    check_eigenpairs(symmetric, both);
}

TEST(LanczosCsr, DiagonalOperatorsAreExact) {
    for (bool explicit_zeros : {false, true}) {
        constexpr int n = 64, k = 3;
        Csr<C> zero;
        for (int row = 0; row < n; ++row) {
            if (explicit_zeros) { zero.add(row, C{1.}); zero.add(row, C{-1.}); }
            zero.end_row();
        }
        for (std::string const which : {"SA", "LA", "LM"})
            for (bool vectors : {false, true}) {
                auto const r = alps::ed::lanczos_csr(zero.offsets, zero.indices, zero.values, n, k, which, 1e-12, 1000, vectors);
                EXPECT_EQ(r.nconv, k);
                EXPECT_EQ(r.eigenvalues, std::vector<double>(k, 0.));
                ASSERT_EQ(r.eigenvectors.size(), vectors ? n * k : 0);
                if (vectors) {
                    EXPECT_EQ(r.residual_norms, std::vector<double>(k, 0.));
                    for (int col = 0; col < k; ++col)
                        for (int row = 0; row < n; ++row)
                            EXPECT_EQ(r.eigenvectors[static_cast<std::size_t>(col * n + row)], C{row == col ? 1. : 0.});
                }
            }
    }
    Csr<double> diagonal;
    for (double value : {-9., -3., -.2, .7, 2., 6.}) { diagonal.add(diagonal.dim, value); diagonal.end_row(); }
    for (std::string const which : {"SA", "LA", "LM"}) {
        auto const r = alps::ed::lanczos_csr(diagonal.offsets, diagonal.indices, diagonal.values, 6, 2, which, 1e-12, 100, true);
        std::vector<double> const expected = which == "SA" ? std::vector<double>{-9., -3.}
            : which == "LA" ? std::vector<double>{2., 6.} : std::vector<double>{-9., 6.};
        EXPECT_EQ(r.eigenvalues, expected) << which;
        check_eigenpairs(diagonal, r, 1e-14);
    }
}

TEST(LanczosCsr, DegenerateMultipletsAreResolvedByDeflation) {
    // Two identical random tridiagonal blocks: every eigenvalue is doubly
    // degenerate. The lowest five must come back in pairs and agree with
    // the dense spectrum of the full matrix.
    constexpr int block = 150, n = 2 * block, k = 5;
    std::mt19937_64 engine(7);
    std::uniform_real_distribution<double> random(-1., 1.);
    std::vector<double> diagonal(block), hopping(block - 1);
    for (auto& d : diagonal) d = random(engine);
    for (auto& t : hopping) t = random(engine);
    Csr<double> m;
    std::vector<double> dense(static_cast<std::size_t>(n) * n, 0.);
    for (int b = 0; b < 2; ++b)
        for (int i = 0; i < block; ++i) {
            int const row = b * block + i;
            if (i > 0) { m.add(row - 1, hopping[static_cast<std::size_t>(i - 1)]); dense[static_cast<std::size_t>((row - 1) * n + row)] = hopping[static_cast<std::size_t>(i - 1)]; }
            m.add(row, diagonal[static_cast<std::size_t>(i)]); dense[static_cast<std::size_t>(row * n + row)] = diagonal[static_cast<std::size_t>(i)];
            if (i + 1 < block) { m.add(row + 1, hopping[static_cast<std::size_t>(i)]); dense[static_cast<std::size_t>((row + 1) * n + row)] = hopping[static_cast<std::size_t>(i)]; }
            m.end_row();
        }
    auto const reference = alps::ed::dense_eigensystem(dense, n, false, k, "SA").eigenvalues;
    auto const r = alps::ed::lanczos_csr(m.offsets, m.indices, m.values, n, k);
    ASSERT_EQ(r.eigenvalues.size(), k);
    for (int i = 0; i < k; ++i) EXPECT_NEAR(r.eigenvalues[static_cast<std::size_t>(i)], reference[static_cast<std::size_t>(i)], 1e-8) << i;
    EXPECT_NEAR(r.eigenvalues[0], r.eigenvalues[1], 1e-8);
    EXPECT_NEAR(r.eigenvalues[2], r.eigenvalues[3], 1e-8);
    EXPECT_GT(r.runs, 2);
    check_eigenpairs(m, r);
    auto const top = alps::ed::lanczos_csr(m.offsets, m.indices, m.values, n, 2, "LA");
    auto const highest = alps::ed::dense_eigensystem(dense, n, false, 2, "LA").eigenvalues;
    ASSERT_EQ(top.eigenvalues.size(), 2);
    EXPECT_NEAR(top.eigenvalues[0], highest[0], 1e-8);
    EXPECT_NEAR(top.eigenvalues[1], highest[1], 1e-8);
    check_eigenpairs(m, top);
}

TEST(LanczosCsr, DuplicateAndReversedEntriesDescribeTheSameOperator) {
    auto const m = chain<double>(120, -1.);
    Csr<double> split;
    for (std::int64_t row = 0; row < m.dim; ++row) {
        for (auto pos = m.offsets[static_cast<std::size_t>(row) + 1]; pos-- > m.offsets[static_cast<std::size_t>(row)];) {
            split.add(m.indices[static_cast<std::size_t>(pos)], .25 * m.values[static_cast<std::size_t>(pos)]);
            split.add(m.indices[static_cast<std::size_t>(pos)], .75 * m.values[static_cast<std::size_t>(pos)]);
        }
        split.end_row();
    }
    auto const original = alps::ed::lanczos_csr(m.offsets, m.indices, m.values, m.dim, 2);
    auto const duplicated = alps::ed::lanczos_csr(split.offsets, split.indices, split.values, split.dim, 2);
    ASSERT_EQ(original.eigenvalues.size(), 2);
    ASSERT_EQ(duplicated.eigenvalues.size(), 2);
    for (std::size_t i = 0; i < 2; ++i) EXPECT_NEAR(original.eigenvalues[i], duplicated.eigenvalues[i], 1e-9);
    check_eigenpairs(split, duplicated);
}

TEST(LanczosCsr, ContractAndValidation) {
    auto const m = chain<double>(40, -1.);
    auto solve = [&](auto const& offsets, auto const& indices, auto const& values, std::int64_t dim, int k,
                     std::string const& which = "SA", double tol = 0.) {
        return alps::ed::lanczos_csr(offsets, indices, values, dim, k, which, tol, 0, false);
    };
    EXPECT_THROW(solve(m.offsets, m.indices, m.values, m.dim, 40), std::invalid_argument);
    EXPECT_THROW(solve(m.offsets, m.indices, m.values, m.dim, 0), std::invalid_argument);
    EXPECT_THROW(solve(m.offsets, m.indices, m.values, m.dim, 1, "XX"), std::invalid_argument);
    EXPECT_THROW(solve(m.offsets, m.indices, m.values, m.dim, 1, "SM"), std::invalid_argument);
    EXPECT_THROW(solve(std::vector<std::int64_t>{}, m.indices, m.values, m.dim, 1), std::invalid_argument);
    EXPECT_THROW(solve(m.offsets, std::vector<std::int64_t>(m.indices.begin(), m.indices.end() - 1), m.values, m.dim, 1),
                 std::invalid_argument);
    auto bad_offsets = m.offsets;
    bad_offsets[0] = 1;
    EXPECT_THROW(solve(bad_offsets, m.indices, m.values, m.dim, 1), std::invalid_argument);
    bad_offsets = m.offsets;
    bad_offsets.back() += 1;
    EXPECT_THROW(solve(bad_offsets, m.indices, m.values, m.dim, 1), std::invalid_argument);
    bad_offsets = m.offsets;
    bad_offsets[2] = -1;
    EXPECT_THROW(solve(bad_offsets, m.indices, m.values, m.dim, 1), std::invalid_argument);
    auto bad_indices = m.indices;
    bad_indices[0] = m.dim;
    EXPECT_THROW(solve(m.offsets, bad_indices, m.values, m.dim, 1), std::invalid_argument);
    bad_indices[0] = -1;
    EXPECT_THROW(solve(m.offsets, bad_indices, m.values, m.dim, 1), std::invalid_argument);
    auto bad_values = m.values;
    bad_values[0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_THROW(solve(m.offsets, m.indices, bad_values, m.dim, 1), std::invalid_argument);
    EXPECT_THROW(solve(m.offsets, m.indices, m.values, m.dim, 1, "SA", std::numeric_limits<double>::infinity()),
                 std::invalid_argument);
    EXPECT_THROW(solve(std::vector<std::int64_t>{}, std::vector<std::int64_t>{}, std::vector<double>{},
                       std::int64_t{std::numeric_limits<int>::max()} + 1, 1), std::invalid_argument);
}

TEST(LanczosCsr, NonConvergenceIsReported) {
    // Two Lanczos steps cannot converge the lowest level of a 500-site chain
    // to 1e-14. Failure must be explicit, never a partial result.
    auto const m = chain<double>(500, -1., 0.3);
    try {
        (void)alps::ed::lanczos_csr(m.offsets, m.indices, m.values, m.dim, 3, "SA", 1e-14, 2, true);
        FAIL() << "an inadequate step limit must report non-convergence";
    } catch (std::runtime_error const& error) {
        EXPECT_NE(std::string(error.what()).find("did not converge"), std::string::npos) << error.what();
    }
}
