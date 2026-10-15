// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/tests/models/sparse_from_local_lock.cpp.
// Adaptation: <alps/ed/...> include names; converted to GoogleTest; spectra
// from alps::ed::dense_eigensystem and the models test helper instead of Eigen.
//
// Lock for alps::ed::build_sparse_from_local: the CSR assembled from a kit
// descriptor must be spectrum-identical to the kit's dense reference assembly
// (src/alps/models/tests/dense_e0.hpp — an independent Kronecker-product
// implementation with the opposite site significance, so index conventions
// cannot silently cancel), across spin / fermion-with-JW / mixed-dimension
// Kondo models; QN-sector constraints must project to the right dimensions
// and preserve the ground state of conserving Hamiltonians; and the exact
// 2-site Heisenberg anchor pins the absolute scale.

#include "dense_e0.hpp"

#include <alps/ed/eigensolver.hpp>
#include <alps/ed/sparse_from_local.hpp>
#include <alps/models/builtin.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace lat = alps::lattice;
namespace mdl = alps::models;

namespace {

std::map<int, std::vector<double>> parities(mdl::Hamiltonian const& H) {
    std::map<int, std::vector<double>> out;
    for (int i = 0; i < H.num_sites(); ++i) {
        int const t = H.lat().sites()[static_cast<std::size_t>(i)].site_type;
        if (!out.count(t)) out[t] = H.site_for(i)->parity_operator;
    }
    return out;
}

std::vector<double> csr_spectrum(alps::ed::SparseFromLocal const& m) {
    auto const n = static_cast<std::size_t>(m.dim);
    std::vector<double> dense(n * n, 0.);   // column-major
    for (std::int64_t r = 0; r < m.dim; ++r)
        for (std::int64_t p = m.indptr[static_cast<std::size_t>(r)];
             p < m.indptr[static_cast<std::size_t>(r) + 1]; ++p)
            dense[static_cast<std::size_t>(m.indices[static_cast<std::size_t>(p)]) * n
                  + static_cast<std::size_t>(r)] = m.data[static_cast<std::size_t>(p)];
    // The kernel promises a symmetric real matrix for these models.
    double asym = 0.;
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t c = 0; c < n; ++c)
            asym = std::max(asym, std::abs(dense[c * n + r] - dense[r * n + c]));
    EXPECT_LE(asym, 1e-12) << "CSR not symmetric";
    return alps::ed::dense_eigensystem(std::move(dense), n, false).eigenvalues;
}

std::vector<double> dense_spectrum(mdl::Hamiltonian const& H,
                                   alps::models::CompiledModel const& d) {
    return alps_test::hermitian_eigenvalues(alps_test::dense_hamiltonian(H, d));
}

void check_full_spectrum(std::string const& name, mdl::Hamiltonian const& H) {
    SCOPED_TRACE(name);
    auto const d = H.to_local_decomposition();
    auto const csr = alps::ed::build_sparse_from_local(d, parities(H));
    auto const a = csr_spectrum(csr);
    auto const b = dense_spectrum(H, d);
    ASSERT_EQ(a.size(), b.size());
    double const scale =
        std::max(1.0, std::max(std::abs(a.front()), std::abs(a.back())));
    for (std::size_t k = 0; k < a.size(); ++k)
        ASSERT_NEAR(a[k], b[k], 1e-10 * scale) << "eigenvalue " << k;
}

lat::Lattice kondo_geometry(int L) {
    lat::UnitCell uc(1);
    uc.add_vertex({0.0}, 0);
    uc.add_vertex({0.5}, 1);
    uc.add_edge(0, 0, 0, {}, {1});
    uc.add_edge(0, 1, 1);
    return uc.tile({L});
}

} // namespace

TEST(SparseFromLocalLock, CooTermTraversalContract) {
    // Pin the unreduced contribution stream across scalar, single-output,
    // zero-column and branching factors. Projection and JW zero weights
    // must not reorder the contributions that remain.
    alps::models::CompiledOperator d;
    d.site_dims = {2, 2};
    d.terms = {{{0., 2.}, {}, {}},
        {{1., .5}, {{0, {2., 0., 0., -1.}}, {1, {0., 3., 4., 0.}}}, {0}},
        {{-1., -.5}, {{0, {2., 0., 0., -1.}}, {1, {0., 3., 4., 0.}}}, {0}},
        {1., {{0, {0., 0., 5., 0.}}}, {}},
        {.25, {{0, {1., 2., 3., 4.}}}, {1}}};
    std::vector<std::int64_t> const states{1, 0, 0, 1, 1, 1};
    std::vector<std::vector<double>> const parity{{1., -1.}, {0., -1.}};
    auto const matrix = alps::ed::build_coo_from_local(d, states, 3, parity);
    EXPECT_EQ(matrix.rows, (std::vector<std::int64_t>{0, 1, 2, 2, 0, 2, 0, 2, 1, 2, 1, 2}));
    EXPECT_EQ(matrix.cols, (std::vector<std::int64_t>{0, 1, 2, 0, 2, 0, 2, 1, 1, 1, 2, 2}));
    EXPECT_EQ(matrix.data, (std::vector<std::complex<double>>{
                  {0., 2.}, {0., 2.}, {0., 2.}, {4., 2.}, {3., 1.5},
                  {-4., -2.}, {-3., -1.5}, 5., -.25, -.75, -.5, -1.}));
    EXPECT_TRUE(matrix.is_complex);
}

TEST(SparseFromLocalLock, ComplexExplicitBasis) {
    alps::models::CompiledOperator d;
    d.site_dims = {2};
    d.terms = {{1., {{0, {0., {0., 1e-15}, 2., 0.}}}, {}}};
    std::vector<std::int64_t> const states{1, 0};
    auto const matrix = alps::ed::build_coo_from_local(d, states, 2);
    EXPECT_EQ(matrix.rows, (std::vector<std::int64_t>{1, 0}));
    EXPECT_EQ(matrix.cols, (std::vector<std::int64_t>{0, 1}));
    EXPECT_EQ(matrix.data, (std::vector<std::complex<double>>{{0., 1e-15}, 2.}));
    EXPECT_TRUE(matrix.is_complex);
    EXPECT_TRUE(alps::ed::build_coo_from_local(d, {}, 0).data.empty());
    auto rejected = [&](auto const& descriptor, auto const& basis, std::size_t count) {
        try { alps::ed::build_coo_from_local(descriptor, basis, count); }
        catch (std::invalid_argument const&) { return true; }
        return false;
    };
    EXPECT_TRUE(rejected(d, states, 3)) << "invalid basis shape accepted";
    EXPECT_TRUE(rejected(d, std::vector<std::int64_t>{0, 0}, 2)) << "duplicate states accepted";
    EXPECT_TRUE(rejected(d, std::vector<std::int64_t>{-1}, 1)) << "negative local state accepted";
    EXPECT_TRUE(rejected(d, std::vector<std::int64_t>{2}, 1)) << "out-of-range local state accepted";
    d.site_dims[0] = 0;
    EXPECT_TRUE(rejected(d, states, 2)) << "zero local dimension accepted";
    d.site_dims[0] = 2;
    d.terms[0].factors[0].site = 1;
    EXPECT_TRUE(rejected(d, states, 2)) << "out-of-range factor site accepted";
    d.terms[0].factors[0].site = 0;
    d.terms[0].jordan_wigner = {1};
    EXPECT_TRUE(rejected(d, states, 2)) << "out-of-range JW site accepted";
    d.terms[0].jordan_wigner.clear();
    d.terms[0].factors[0].matrix.pop_back();
    EXPECT_TRUE(rejected(d, states, 2)) << "short operator accepted";
}

TEST(SparseFromLocalLock, UnencodableExplicitBasis) {
    alps::models::CompiledOperator d;
    d.site_dims.assign(70, 2);
    d.terms = {{{-.25, -.5},
        {{0, {0., 1., 0., 0.}}, {69, {0., 0., 1., 0.}}}, {}}};
    for (int site = 0; site < 69; ++site) d.terms[0].jordan_wigner.push_back(site);
    std::vector<std::int64_t> states(140, 0);
    states[0] = states[35] = states[70 + 35] = states[70 + 69] = 1;
    std::vector<std::vector<double>> const parity(70, {1., -1.});
    auto const matrix = alps::ed::build_coo_from_local(d, states, 2, parity);
    EXPECT_EQ(matrix.rows, std::vector<std::int64_t>{1});
    EXPECT_EQ(matrix.cols, std::vector<std::int64_t>{0});
    EXPECT_EQ(matrix.data, (std::vector<std::complex<double>>{{-.25, -.5}}));
    // 3 * 2^62 states fit unsigned codes, but this selected input code
    // exceeds INT64_MAX. Replacing its high digit must not overflow.
    d.site_dims.assign(63, 2);
    d.site_dims[0] = 3;
    d.terms = {{1., {{0, {0., 0., 1., 0., 0., 0., 0., 0., 0.}}}, {}}};
    states.assign(126, 0);
    states[0] = 2;
    auto const high_code = alps::ed::build_coo_from_local(d, states, 2);
    EXPECT_EQ(high_code.rows, std::vector<std::int64_t>{1});
    EXPECT_EQ(high_code.cols, std::vector<std::int64_t>{0});
    EXPECT_EQ(high_code.data, std::vector<std::complex<double>>{1.});
    d.site_dims.clear();
    d.terms = {{{0., 2.}, {}, {}}};
    EXPECT_EQ(alps::ed::build_coo_from_local(d, {}, 1).data,
              (std::vector<std::complex<double>>{{0., 2.}}));
}

// Matrix entries, rather than spectra alone, pin the native caller's
// input-row convention, site significance, duplicate accumulation and
// explicit zeros.
TEST(SparseFromLocalLock, DescriptorMatrixContract) {
    alps::models::CompiledModel d;
    d.op.site_dims = {2, 3};
    d.site_types = {0, 1};
    d.op.terms = {{2., {{1, {0., 1., 0., 0., 0., 3., 4., 0., 0.}}}, {0}},
        {-.5, {}, {}},
        {.5, {}, {}}};
    auto const m = alps::ed::build_sparse_from_local(d, {{0, {1., -1.}}});
    EXPECT_EQ(m.states, (std::vector<std::uint64_t>{0, 1, 2, 3, 4, 5}));
    EXPECT_EQ(m.indptr, (std::vector<std::int64_t>{0, 2, 4, 6, 8, 10, 12}));
    EXPECT_EQ(m.indices, (std::vector<std::int64_t>{0, 4, 1, 5, 0, 2, 1, 3, 2, 4, 3, 5}));
    EXPECT_EQ(m.data, (std::vector<double>{0., 8., 0., -8., 2., 0., -2., 0., 6., 0., -6., 0.}));
    d.op.terms[0].coefficient = {2., 1e-15};
    EXPECT_EQ(alps::ed::build_sparse_from_local(d).data.size(), 12);
    d.op.terms[0].coefficient = {2., 1e-3};
    EXPECT_THROW(alps::ed::build_sparse_from_local(d), std::invalid_argument);
}

TEST(SparseFromLocalLock, HeisenbergChain2Anchor) {
    auto const H = mdl::heisenberg(lat::chain(2, 1.0, lat::Boundary::open));
    auto const csr = alps::ed::build_sparse_from_local(
        H.to_local_decomposition(), parities(H));
    ASSERT_EQ(csr.dim, 4);
    auto const ev = csr_spectrum(csr);
    EXPECT_NEAR(ev[0], -0.75, 1e-12);
    for (std::size_t k = 1; k < ev.size(); ++k)
        EXPECT_NEAR(ev[k], 0.25, 1e-12) << "E" << k;
}

TEST(SparseFromLocalLock, SpectrumParityWithTheDenseAssembly) {
    check_full_spectrum("heisenberg_chain4", mdl::heisenberg(lat::chain(4)));
    check_full_spectrum("xxz_ladder3", mdl::xxz(lat::ladder(3), 0.6, 1.0));
    check_full_spectrum("tfi_chain4", mdl::transverse_field_ising(lat::chain(4), 1.0, 0.7));
    check_full_spectrum("hubbard_chain3", mdl::hubbard(lat::chain(3), 1.0, 2.0, 0.3));
    check_full_spectrum("spinless_chain4", mdl::spinless_fermion_model(lat::chain(4), 1.0, 0.2, 0.5));
    check_full_spectrum("kondo_L2", mdl::kondo_lattice(kondo_geometry(2), 1.0, 1.3, 0.2));
}

TEST(SparseFromLocalLock, QuantumNumberSectors) {
    {
        auto const H = mdl::heisenberg(lat::chain(4));
        auto const d = H.to_local_decomposition();
        auto const full = alps::ed::build_sparse_from_local(d, parities(H));
        auto const sz0 = alps::ed::build_sparse_from_local(d, parities(H), {{"Sz", 0.0}});
        auto const sz1 = alps::ed::build_sparse_from_local(d, parities(H), {{"Sz", 1.0}});
        EXPECT_EQ(sz0.dim, 6);
        EXPECT_EQ(sz1.dim, 4);
        EXPECT_NEAR(csr_spectrum(full).front(), csr_spectrum(sz0).front(), 1e-12);
    }
    {
        // Hubbard half filling: N = 2, Sz = 0 on a 2-site chain has exactly
        // the 4 states {(ud,0),(0,ud),(u,d),(d,u)} and holds the GS.
        auto const H = mdl::hubbard(lat::chain(2), 1.0, 4.0, 0.0);
        auto const d = H.to_local_decomposition();
        auto const sec = alps::ed::build_sparse_from_local(d, parities(H), {{"N", 2.0}, {"Sz", 0.0}});
        EXPECT_EQ(sec.dim, 4);
        auto const full = alps::ed::build_sparse_from_local(d, parities(H));
        EXPECT_NEAR(csr_spectrum(full).front(), csr_spectrum(sec).front(), 1e-12);
    }
}
