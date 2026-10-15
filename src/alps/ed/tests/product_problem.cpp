// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/tests/product_problem.cpp.
// Adaptation: C++17 (no designated initializers or std::numbers), <alps/ed/...>
// include names; the Lanczos cases use larger sectors and the IETL driver's
// controls; expectation values added.
#include <alps/ed/product_problem.hpp>
#include <alps/models/builtin.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <future>
#include <limits>
#include <stdexcept>

namespace {
using namespace alps::ed;
namespace mdl = alps::models;
using C = std::complex<double>;
double const pi = std::acos(-1.);
auto chain(int n) { return alps::lattice::chain(n, 1., alps::lattice::Boundary::open); }
EigenstateOptions options(int count, EigenstateMethod method = EigenstateMethod::automatic,
                          std::string which = "SA") {
    EigenstateOptions o;
    o.count = count;
    o.method = method;
    o.which = std::move(which);
    return o;
}
ProductOptions sector(std::map<std::string, double> charges) {
    ProductOptions o;
    o.sector = std::move(charges);
    return o;
}
std::vector<C> matrix(ProductProblem const& problem) {
    auto n = problem.dimension();
    std::vector<C> result(n*n);
    std::visit([&](auto const& values) {
        for (std::size_t row = 0; row < n; ++row)
            for (auto entry = problem.row_offsets()[row]; entry < problem.row_offsets()[row+1]; ++entry)
                result[row*n+static_cast<std::size_t>(problem.column_indices()[static_cast<std::size_t>(entry)])]
                    = values[static_cast<std::size_t>(entry)];
    }, problem.matrix_values());
    return result;
}
void residual(ProductEigenResult const& result, double tolerance = 1e-10) {
    auto const n = result.problem.dimension(), k = result.eigenvalues.size();
    auto const h = matrix(result.problem);
    std::visit([&](auto const& vectors) {
        ASSERT_EQ(vectors.size(), n*k);
        for (std::size_t col = 0; col < k; ++col) {
            double error = 0.;
            for (std::size_t row = 0; row < n; ++row) {
                C value = -result.eigenvalues[col]*vectors[col*n+row];
                for (std::size_t j = 0; j < n; ++j) value += h[row*n+j]*vectors[col*n+j];
                error += std::norm(value);
            }
            EXPECT_LT(std::sqrt(error), tolerance);
            if (result.method == EigenstateMethod::lanczos)
                EXPECT_NEAR(result.residual_norms[col], std::sqrt(error), 1e-12);
            for (std::size_t other = 0; other < k; ++other) {
                C overlap = 0.;
                for (std::size_t row = 0; row < n; ++row)
                    overlap += std::conj(C(vectors[col*n+row]))*vectors[other*n+row];
                EXPECT_NEAR(std::abs(overlap - C(col == other)), 0., 1e-10);
            }
        }
    }, result.eigenvectors);
}
}
TEST(ProductProblem, NativeModelWorkflowAndDetachedResult) {
    auto result = [] {
        auto model = mdl::heisenberg(chain(2));
        auto prepared = prepare_product(model, sector({{"Sz",0.}}));
        model.add(100.*mdl::TermSum::located(mdl::spin(), "Sz", 0));
        EXPECT_EQ(prepared.dimension(), 2);
        EXPECT_EQ(prepared.states(), (std::vector<std::int32_t>{0,1,1,0}));
        EXPECT_EQ(matrix(prepared), (std::vector<C>{-.25,.5,.5,-.25}));
        return prepared.eigenstates(options(8));
    }();
    ASSERT_EQ(result.eigenvalues.size(), 2);
    EXPECT_NEAR(result.eigenvalues[0], -.75, 1e-14);
    EXPECT_NEAR(result.eigenvalues[1], .25, 1e-14);
    EXPECT_TRUE(std::holds_alternative<std::vector<double>>(result.eigenvectors));
    EXPECT_EQ(result.method, EigenstateMethod::dense);
    EXPECT_FALSE(result.tolerance);
    EXPECT_EQ(result.problem.sector(), (std::map<std::string, double>{{"Sz", 0.}}));
    ASSERT_EQ(result.problem.physical_basis().sites.size(), 2);
    EXPECT_EQ(result.problem.physical_basis().conserved_quantum_numbers, (std::vector<std::string>{"Sz"}));
    residual(result);
    auto original = result.problem;
    auto moved = std::move(original);
    EXPECT_THROW(original.eigenstates(), std::logic_error);
    EXPECT_EQ(moved.eigenstates().eigenvalues, result.eigenvalues);
}
TEST(ProductProblem, ValidationPrecedesProjectionAndUnusedChargesStayOptional) {
    auto spin = mdl::spin();
    auto flip = mdl::TermSum::located(spin, "Sx", 0);
    EXPECT_NO_THROW(prepare_product(flip, {spin}, {"Sz"}));
    EXPECT_THROW(prepare_product(flip, {spin}, {"Sz"}, mdl::FermionConvention::strict,
                                sector({{"Sz",99.}})), std::invalid_argument);
    auto raising = mdl::TermSum::located(spin, "Sp", 0);
    EXPECT_THROW(prepare_product(raising, {spin}), std::invalid_argument);
    EXPECT_THROW(prepare_product(raising, {spin}, {"Sz"}, mdl::FermionConvention::strict,
                                sector({{"Sz",99.}})), std::invalid_argument);
    EXPECT_THROW(prepare_product(flip, {spin}, {}, mdl::FermionConvention::strict,
                                sector({{"Sz",0.}})), std::invalid_argument);
    auto empty = prepare_product(mdl::heisenberg(chain(2)), sector({{"Sz",99.}}));
    EXPECT_EQ(empty.dimension(), 0);
    auto result = empty.eigenstates(options(1, EigenstateMethod::lanczos));
    EXPECT_TRUE(result.eigenvalues.empty());
    EXPECT_EQ(result.method, EigenstateMethod::dense);
    EXPECT_TRUE(result.sparse_requested);
}
TEST(ProductProblem, ComplexFermionSelectionAndResiduals) {
    auto site = mdl::spinless_fermion();
    constexpr int L = 300;
    mdl::Hamiltonian model(chain(L), site, mdl::FermionConvention::strict, {"N"});
    model.add_bond(-std::polar(1., .3), [site](int i, int j) {
        return mdl::TermSum::located(site, "c_dag", i)*mdl::TermSum::located(site, "c", j);
    }, {}, true);
    model.add_site(.125, [site](int i) { return mdl::TermSum::located(site, "n", i); });
    auto prepared = prepare_product(model, sector({{"N",1.}}));
    ASSERT_EQ(prepared.dimension(), L);
    ASSERT_TRUE(std::holds_alternative<std::vector<C>>(prepared.matrix_values()));
    for (auto which : {"SA", "LA", "LM"}) {
        SCOPED_TRACE(which);
        auto dense = prepared.eigenstates(options(3, EigenstateMethod::dense, which));
        auto sparse = prepared.eigenstates(options(3, EigenstateMethod::lanczos, which));
        ASSERT_EQ(sparse.eigenvalues.size(), 3);
        ASSERT_EQ(dense.eigenvalues.size(), 3);
        for (std::size_t i = 0; i < 3; ++i) EXPECT_NEAR(dense.eigenvalues[i], sparse.eigenvalues[i], 1e-9);
        EXPECT_EQ(sparse.method, EigenstateMethod::lanczos);
        EXPECT_EQ(sparse.max_iterations, 1000);
        EXPECT_EQ(sparse.tolerance, 1e-10);
        EXPECT_TRUE(sparse.iterations && *sparse.iterations > 0);
        residual(dense); residual(sparse, 1e-7);
    }
    auto dense_sm = prepared.eigenstates(options(2, EigenstateMethod::dense, "SM"));
    auto automatic_sm = prepared.eigenstates(options(2, EigenstateMethod::automatic, "SM"));
    EXPECT_EQ(automatic_sm.method, EigenstateMethod::dense);
    EXPECT_EQ(automatic_sm.eigenvalues, dense_sm.eigenvalues);
    EXPECT_THROW(prepared.eigenstates(options(2, EigenstateMethod::lanczos, "SM")), std::invalid_argument);
    auto ground = prepared.eigenstates(options(1));
    EXPECT_EQ(ground.method, EigenstateMethod::dense);   // automatic: dimension <= 1000
    EXPECT_NEAR(ground.eigenvalues[0], .125-2*std::cos(pi/(L+1)), 1e-12);
}
TEST(ProductProblem, HugeFullSpaceAndNeutralLocalCharges) {
    auto spin = mdl::spin();
    std::vector<mdl::SiteBasisPtr> sites(70, spin);
    auto prepared = prepare_product(mdl::TermSum::located(spin,"Sz",0), sites, {"Sz"},
        mdl::FermionConvention::strict, sector({{"Sz",35.}}));
    ASSERT_EQ(prepared.dimension(),1);
    EXPECT_EQ(prepared.states().size(),70);
    EXPECT_NEAR(prepared.eigenstates().eigenvalues[0],.5,1e-14);
    auto neutral = std::make_shared<mdl::SiteBasis>(*spin);
    neutral->quantum_numbers.erase("Sz");
    auto mixed = prepare_product(mdl::TermSum::located(spin,"Sz",0), {spin, neutral}, {"Sz"},
        mdl::FermionConvention::strict, sector({{"Sz",.5}}));
    EXPECT_EQ(mixed.dimension(),2);
    neutral->basis_labels = {"edited", "caller"};
    EXPECT_NE(mixed.physical_basis().sites[1].labels, neutral->basis_labels);
    EXPECT_EQ(mixed.physical_basis().sites[1].quantum_numbers.at("Sz"), (std::vector<double>{0.,0.}));
}
TEST(ProductProblem, ControlsFailuresAndConcurrentReuse) {
    auto small = prepare_product(mdl::heisenberg(chain(8)), sector({{"Sz",0.}}));
    EXPECT_THROW(small.eigenstates(options(0)), std::invalid_argument);
    EXPECT_THROW(small.eigenstates(options(1, EigenstateMethod::automatic, "invalid")), std::invalid_argument);
    auto zero_tolerance = options(1);
    zero_tolerance.tolerance = 0.;
    EXPECT_THROW(small.eigenstates(zero_tolerance), std::invalid_argument);
    EXPECT_THROW(small.eigenstates(options(1, static_cast<EigenstateMethod>(-1))), std::invalid_argument);
    // Dimension 70 is below the Lanczos threshold: a requested Lanczos solve is dense.
    auto fallback = small.eigenstates(options(3, EigenstateMethod::lanczos));
    EXPECT_EQ(fallback.method, EigenstateMethod::dense);
    EXPECT_TRUE(fallback.sparse_requested);
    residual(fallback);

    auto prepared = prepare_product(mdl::heisenberg(chain(12)), sector({{"Sz",0.}}));
    ASSERT_EQ(prepared.dimension(), 924);
    auto starved = options(3, EigenstateMethod::lanczos);
    starved.max_iterations = 2;
    starved.tolerance = 1e-14;
    EXPECT_THROW(prepared.eigenstates(starved), std::runtime_error);
    auto work = [prepared] { return prepared.eigenstates(options(3, EigenstateMethod::lanczos)); };
    auto first = std::async(std::launch::async, work), second = std::async(std::launch::async, work);
    auto a = first.get(), b = second.get();
    EXPECT_EQ(a.eigenvalues, b.eigenvalues);
    EXPECT_EQ(a.method, EigenstateMethod::lanczos);
    residual(a, 1e-7); residual(b, 1e-7);
    auto dense = prepared.eigenstates(options(3, EigenstateMethod::dense));
    for (std::size_t i = 0; i < 3; ++i) EXPECT_NEAR(a.eigenvalues[i], dense.eigenvalues[i], 1e-8);
}
TEST(ProductProblem, ObservablesOnEigenvectors) {
    auto const spin = mdl::spin();
    auto const model = mdl::heisenberg(chain(6));
    auto const prepared = prepare_product(model, sector({{"Sz",0.}}));
    auto const result = prepared.eigenstates(options(2, EigenstateMethod::dense));
    // Total Sz is fixed by the sector; the nearest-neighbour correlation of
    // the singlet ground state is E0/(L-1) per bond.
    mdl::TermSum total;
    for (int i = 0; i < 6; ++i) total = total + mdl::TermSum::located(spin, "Sz", i);
    auto const sz = expectation_values(result, mdl::compile_operator(total, std::vector<mdl::SiteBasisPtr>(6, spin)));
    ASSERT_EQ(sz.size(), 2);
    for (auto value : sz) EXPECT_NEAR(std::abs(value), 0., 1e-12);
    auto const bond_energy = expectation_values(result, model.operator_decomposition().op);
    ASSERT_EQ(bond_energy.size(), 2);
    EXPECT_NEAR(bond_energy[0].real(), result.eigenvalues[0], 1e-12);
    EXPECT_NEAR(bond_energy[1].real(), result.eigenvalues[1], 1e-12);
    // A charge-changing operator leaves the sector: its projection vanishes.
    auto const raising = expectation_values(result,
        mdl::compile_operator(mdl::TermSum::located(spin, "Sp", 2), std::vector<mdl::SiteBasisPtr>(6, spin)));
    for (auto value : raising) EXPECT_EQ(value, C{});
    alps::models::CompiledOperator wrong;
    wrong.site_dims = {2, 2};
    EXPECT_THROW(prepared.operator_matrix(wrong), std::invalid_argument);
}
