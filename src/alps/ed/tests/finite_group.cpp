// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/tests/finite_group.cpp.
// Adaptation: C++17 (no designated initializers or std::numbers),
// <alps/ed/...> include names; the dense checks use a small local matrix
// helper (alpsx used Eigen).
#include <alps/ed/product_problem.hpp>
#include <alps/ed/projected_problem.hpp>
#include <alps/models/builtin.hpp>
#include <alps/models/physical_basis.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <future>
#include <stdexcept>
#include <vector>

namespace {
using namespace alps::ed;
namespace mdl = alps::models;
using C = std::complex<double>;
double const pi = std::acos(-1.);

// Minimal dense complex matrix for independent checks.
struct Dense {
    std::size_t rows = 0, cols = 0;
    std::vector<C> data; // row-major
    Dense() = default;
    Dense(std::size_t r, std::size_t c) : rows(r), cols(c), data(r * c) {}
    static Dense identity(std::size_t n) {
        Dense m(n, n);
        for (std::size_t i = 0; i < n; ++i) m(i, i) = 1.;
        return m;
    }
    C& operator()(std::size_t r, std::size_t c) { return data[r * cols + c]; }
    C operator()(std::size_t r, std::size_t c) const { return data[r * cols + c]; }
    Dense adjoint() const {
        Dense m(cols, rows);
        for (std::size_t r = 0; r < rows; ++r)
            for (std::size_t c = 0; c < cols; ++c) m(c, r) = std::conj((*this)(r, c));
        return m;
    }
    Dense col(std::size_t c) const {
        Dense v(rows, 1);
        for (std::size_t r = 0; r < rows; ++r) v(r, 0) = (*this)(r, c);
        return v;
    }
    double norm() const {
        double sum = 0.;
        for (auto z : data) sum += std::norm(z);
        return std::sqrt(sum);
    }
};
Dense operator*(Dense const& a, Dense const& b) {
    Dense m(a.rows, b.cols);
    for (std::size_t r = 0; r < a.rows; ++r)
        for (std::size_t k = 0; k < a.cols; ++k)
            for (std::size_t c = 0; c < b.cols; ++c) m(r, c) += a(r, k) * b(k, c);
    return m;
}
Dense operator*(double s, Dense a) {
    for (auto& z : a.data) z *= s;
    return a;
}
Dense operator-(Dense a, Dense const& b) {
    for (std::size_t i = 0; i < a.data.size(); ++i) a.data[i] -= b.data[i];
    return a;
}
Dense& operator+=(Dense& a, Dense const& b) {
    for (std::size_t i = 0; i < a.data.size(); ++i) a.data[i] += b.data[i];
    return a;
}
EigenstateOptions options(int count, EigenstateMethod method = EigenstateMethod::automatic) {
    EigenstateOptions o;
    o.count = count;
    o.method = method;
    return o;
}
std::vector<FiniteGroupIrrep> family(std::size_t n) {
    std::vector<FiniteGroupIrrep> result{dihedral_irrep(n, 0, 1), dihedral_irrep(n, 0, -1)};
    for (std::size_t k = 1; k <= n / 2; ++k) {
        result.push_back(dihedral_irrep(n, k));
        if (n % 2 == 0 && k == n / 2) result.push_back(dihedral_irrep(n, k, -1));
    }
    return result;
}
Dense dense(SparseMatrix<C> const& a) {
    Dense result(a.rows, a.columns);
    for (std::size_t r = 0; r < a.rows; ++r)
        for (auto i = a.offsets[r]; i < a.offsets[r + 1]; ++i)
            result(r, static_cast<std::size_t>(a.indices[static_cast<std::size_t>(i)])) = a.values[static_cast<std::size_t>(i)];
    return result;
}
Dense expansion(ProjectedBasis const& basis) {
    Dense result(basis.parent_dimension(), basis.dimension());
    for (std::size_t c = 0; c < basis.dimension(); ++c) {
        std::vector<C> v(basis.dimension());
        v[c] = 1.;
        auto expanded = expand_projected_state(basis, v);
        for (std::size_t r = 0; r < expanded.amplitudes.size(); ++r) result(r, c) = expanded.amplitudes[r];
    }
    return result;
}
// Independent bit-space implementation of a signed site permutation.
Dense action(std::vector<std::size_t> const& permutation, bool fermion) {
    std::size_t const n = permutation.size(), dim = std::size_t{1} << n;
    Dense result(dim, dim);
    for (std::size_t col = 0; col < dim; ++col) {
        std::size_t target = 0;
        int sign = 1;
        for (std::size_t i = 0; i < n; ++i)
            if ((col >> (n - 1 - i)) & 1) {
                target |= std::size_t{1} << (n - 1 - permutation[i]);
                if (fermion)
                    for (std::size_t j = i + 1; j < n; ++j)
                        if (((col >> (n - 1 - j)) & 1) && permutation[i] > permutation[j]) sign = -sign;
            }
        result(target, col) = sign;
    }
    return result;
}
Dense columns(std::vector<C> const& column_major, std::size_t rows, std::size_t count) {
    Dense result(rows, count);
    for (std::size_t c = 0; c < count; ++c)
        for (std::size_t r = 0; r < rows; ++r) result(r, c) = column_major[c * rows + r];
    return result;
}
} // namespace

TEST(FiniteGroup, RelationsIrreducibilityAndCompleteFamilies) {
    for (std::size_t n : {3u, 4u, 5u, 6u}) {
        auto irreps = family(n);
        EXPECT_NO_THROW(validate_irrep_family(irreps));
        for (auto const& irrep : irreps) EXPECT_EQ(irrep.order(), 2 * n);
        irreps.pop_back();
        EXPECT_THROW(validate_irrep_family(irreps), std::invalid_argument);
    }
    EXPECT_THROW((FiniteGroupIrrep({{1, 0}}, {{C(0., 1.)}}, 1, "bad")), std::invalid_argument);
    EXPECT_THROW((FiniteGroupIrrep({{1, 0}}, {{2.}}, 1, "bad")), std::invalid_argument);
    EXPECT_THROW((FiniteGroupIrrep({{1, 1}}, {{1.}}, 1, "bad")), std::invalid_argument);
    EXPECT_THROW((FiniteGroupIrrep({{1, 0}}, {{1., 0., 0., 1.}}, 2, "reducible")), std::invalid_argument);
    EXPECT_THROW((FiniteGroupIrrep({{1, 2, 0}}, {{1.}}, 1, "limit", 2)), std::invalid_argument);
    auto duplicate = family(4);
    duplicate[1] = duplicate[0];
    EXPECT_THROW(validate_irrep_family(duplicate), std::invalid_argument);
    EXPECT_THROW(dihedral_irrep(4, 1, -1), std::invalid_argument);
}

TEST(FiniteGroup, SquareTwoDimensionalIrrepMatchesIndependentCentralProjector) {
    for (bool fermion : {false, true}) {
        std::vector<SpectrumLocalSpace> spaces(4, {{1., fermion ? -1. : 1.}, {}});
        auto irrep = dihedral_irrep(4, 1);
        auto basis = make_finite_group_basis(spaces, irrep);
        auto q = expansion(basis);
        // D4 E characters: +2 at identity, -2 at the half-turn, zero elsewhere.
        Dense expected = .5 * (Dense::identity(16) - action({2, 3, 0, 1}, fermion));
        EXPECT_LT((q * q.adjoint() - expected).norm(), 1e-12);
        EXPECT_LT((q.adjoint() * q - Dense::identity(q.cols)).norm(), 1e-12);
        EXPECT_EQ(basis.dimension() % 2, 0u);
        ASSERT_NE(basis.finite_group(), nullptr);
        EXPECT_EQ(basis.finite_group()->dimension(), 2u);
        Dense all(16, 16);
        for (auto const& item : family(4)) {
            auto columns = expansion(make_finite_group_basis(spaces, item));
            all += columns * columns.adjoint();
        }
        EXPECT_LT((all - Dense::identity(16)).norm(), 1e-12);
    }
}

TEST(FiniteGroup, ComplexCyclicIrrepsMatchMomentumProjectors) {
    std::vector<SpectrumLocalSpace> spaces(5, {{1., -1.}, {{0., 1.}}});
    std::vector<std::size_t> rotation{1, 2, 3, 4, 0};
    for (int k = 0; k < 5; ++k) {
        auto irrep = FiniteGroupIrrep({rotation}, {{std::polar(1., 2 * pi * k / 5)}}, 1, "k");
        auto a = expansion(make_finite_group_basis(spaces, irrep, std::vector<std::optional<double>>{2.}));
        auto b = expansion(make_projected_basis(spaces, std::vector<SpectrumTranslation>{{rotation, 5}},
                                                std::vector<std::optional<double>>{2.}, static_cast<std::size_t>(k)));
        EXPECT_LT((a * a.adjoint() - b * b.adjoint()).norm(), 1e-12);
    }
}

TEST(FiniteGroup, TetrahedralThreeDimensionalIrrepRetainsMultiplicities) {
    // The standard representation of A4 on the vectors orthogonal to (1,1,1,1).
    Dense q(4, 3);
    double const entries[4][3] = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
    for (std::size_t r = 0; r < 4; ++r)
        for (std::size_t c = 0; c < 3; ++c) q(r, c) = entries[r][c] / 2.;
    std::vector<std::vector<std::size_t>> generators{{1, 2, 0, 3}, {1, 0, 3, 2}};
    std::vector<std::vector<C>> matrices;
    for (auto const& g : generators) {
        Dense permutation(4, 4);
        for (std::size_t i = 0; i < 4; ++i) permutation(g[i], i) = 1.;
        Dense d = q.adjoint() * permutation * q;
        auto& values = matrices.emplace_back();
        for (std::size_t r = 0; r < 3; ++r)
            for (std::size_t c = 0; c < 3; ++c) values.push_back(d(r, c));
    }
    FiniteGroupIrrep irrep(generators, matrices, 3, "T");
    ASSERT_EQ(irrep.order(), 12u);
    std::vector<SpectrumLocalSpace> spaces(4, {{1., 1.}, {}});
    auto columns = expansion(make_finite_group_basis(spaces, irrep));
    Dense expected = .25 * (3. * Dense::identity(16) - action({1, 0, 3, 2}, false) - action({2, 3, 0, 1}, false)
                            - action({3, 2, 1, 0}, false));
    EXPECT_LT((columns * columns.adjoint() - expected).norm(), 1e-12);
    EXPECT_EQ(columns.cols % 3, 0u);
}

TEST(FiniteGroup, CompleteSpectraPreserveAllPhysicalMultipletStates) {
    for (int n : {3, 4, 6})
        for (bool fermion : {false, true}) {
            auto model = fermion ? mdl::spinless_fermion_model(alps::lattice::chain(n), .7, .13, .2)
                                 : mdl::heisenberg(alps::lattice::chain(n));
            auto spectrum = finite_spectrum(model, family(static_cast<std::size_t>(n)));
            auto parent = prepare_product(model).eigenstates(options(1000));
            std::vector<double> values;
            for (auto const& block : spectrum.blocks) {
                EXPECT_EQ(block.block_dimension % block.irrep_dimension, 0u);
                EXPECT_EQ(block.group_order, 2 * static_cast<std::size_t>(n));
                values.insert(values.end(), block.eigenvalues.begin(), block.eigenvalues.end());
            }
            std::sort(values.begin(), values.end());
            ASSERT_EQ(values.size(), parent.eigenvalues.size());
            for (std::size_t i = 0; i < values.size(); ++i) EXPECT_NEAR(values[i], parent.eigenvalues[i], 2e-11);
        }
}

TEST(FiniteGroup, DenseLanczosAndCrossIrrepObservablesUsePhysicalProjection) {
    auto model = mdl::heisenberg(alps::lattice::chain(6));
    auto op = mdl::compile_hamiltonian(model).op;
    auto source = prepare_finite(model, dihedral_irrep(6, 1));
    auto target = prepare_finite(model, dihedral_irrep(6, 0));
    auto parent = prepare_projected(model);
    auto q = expansion(source.basis()), r = expansion(target.basis());
    Dense expected = q.adjoint() * dense(parent.hamiltonian()) * q;
    EXPECT_LT((dense(source.hamiltonian()) - expected).norm(), 2e-12);
    EXPECT_EQ(source.basis().finite_group()->label(), "E1");
    for (auto method : {EigenstateMethod::dense, EigenstateMethod::lanczos}) {
        auto result = source.eigenstates(options(2, method));
        EXPECT_EQ(result.spectrum.irrep, "E1");
        EXPECT_EQ(result.spectrum.irrep_dimension, 2u);
        EXPECT_EQ(result.spectrum.group_order, 12u);
        auto v = columns(result.eigenvectors_column_major, source.dimension(), 2);
        for (std::size_t i = 0; i < 2; ++i)
            EXPECT_LT((expected * v.col(i) - result.spectrum.eigenvalues[i] * v.col(i)).norm(), 2e-11);
    }
    mdl::CompiledOperator local{op.site_dims, {{1., {{0, {.5, 0., 0., -.5}}}, {}}}};
    auto actual = dense(projected_sparse_operator(local, source.basis(), target.basis()));
    auto full = dense(projected_sparse_operator(local, parent.basis(), parent.basis()));
    EXPECT_LT((actual - r.adjoint() * full * q).norm(), 2e-12);
    std::vector<C> state(source.dimension(), 0.);
    state[0] = 1.;
    auto applied = apply_projected_operator(local, source.basis(), target.basis(), state);
    EXPECT_LT((columns(applied, applied.size(), 1) - actual.col(0)).norm(), 2e-12);
    auto first = std::async(std::launch::async, [&] { return source.eigenstates(options(2)).spectrum.eigenvalues; });
    EXPECT_EQ(first.get(), source.eigenstates(options(2)).spectrum.eigenvalues);
}

TEST(FiniteGroup, RejectsBrokenSymmetryAndValidatesEmptySectorGlobally) {
    auto model = mdl::heisenberg(alps::lattice::chain(4));
    auto irrep = dihedral_irrep(4, 1);
    EXPECT_EQ(prepare_finite(model, irrep, {{"Sz", 9.}}).dimension(), 0u);
    model.add(.2 * mdl::TermSum::located(mdl::spin(), "Sz", 0));
    EXPECT_THROW(prepare_finite(model, irrep), std::invalid_argument);
    model.add(mdl::TermSum::located(mdl::spin(), "Sp", 0));
    EXPECT_THROW(prepare_finite(model, irrep, {{"Sz", 9.}}), std::invalid_argument);
    auto spaces = std::vector<SpectrumLocalSpace>(4, {{1., 1.}, {}});
    spaces[1].parity = {1., -1.};
    EXPECT_THROW(make_finite_group_basis(spaces, irrep), std::invalid_argument);
}

TEST(FiniteGroup, ScalarOffsetsCannotHideBrokenGroupSymmetry) {
    auto model = mdl::heisenberg(alps::lattice::chain(4));
    auto op = mdl::compile_hamiltonian(model).op;
    auto spaces = projected_spaces(mdl::physical_basis(model));
    op.terms.push_back({1e16, {}, {}});
    op.terms.push_back({.125, {{0, {1., 0., 0., -1.}}}, {}});
    EXPECT_THROW(prepare_finite(op, spaces, dihedral_irrep(4, 1)), std::invalid_argument);
    // Same tensor key, cancellation below a large coefficient's rounding:
    // rounded symbolic aggregation cannot certify this as the zero operator.
    op.terms.clear();
    for (double c : {1e16, .125, -1e16}) op.terms.push_back({c, {{0, {1., 0., 0., -1.}}}, {}});
    EXPECT_THROW(prepare_finite(op, spaces, dihedral_irrep(4, 1)), std::invalid_argument);
    op.terms = {{1e-200, {{0, {1e-200, 0., 0., -1e-200}}}, {}}};
    EXPECT_THROW(prepare_finite(op, spaces, dihedral_irrep(4, 1)), std::invalid_argument);
}

TEST(FiniteGroup, PolarizedHundredSiteSectorDoesNotBuildTheFullHilbertSpace) {
    auto model = mdl::heisenberg(alps::lattice::chain(100));
    auto problem = prepare_finite(model, dihedral_irrep(100, 0), {{"Sz", 50.}});
    ASSERT_EQ(problem.dimension(), 1u);
    EXPECT_DOUBLE_EQ(problem.eigenstates().spectrum.eigenvalues[0], 25.);
    EXPECT_EQ(prepare_finite(model, dihedral_irrep(100, 1), {{"Sz", 50.}}).dimension(), 0u);
}
