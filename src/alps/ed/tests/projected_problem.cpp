// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/tests/projected_problem.cpp.
// Adaptation: C++17 (no designated initializers, abbreviated templates,
// std::span, std::ranges or std::numbers), <alps/ed/...> include names; the
// dispatch and failed-solve cases follow the IETL driver's policy and
// defaults (Lanczos above dimension 1000 and 8*count + 200; tolerance 1e-10,
// min(10*dimension, 1000) steps).
#include <alps/ed/product_problem.hpp>
#include <alps/ed/projected_problem.hpp>
#include <alps/models/builtin.hpp>
#include <alps/models/physical_basis.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <future>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
using namespace alps::ed;
namespace mdl = alps::models;
using C = std::complex<double>;
double const pi = std::acos(-1.);

std::vector<double> sorted(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    return values;
}
template<class F> std::string message(F&& call) {
    try {
        call();
    } catch (std::invalid_argument const& error) {
        return error.what();
    }
    return "no error";
}
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
ProjectedOptions projected(std::map<std::string, double> sector, std::vector<double> momentum = {}) {
    ProjectedOptions o;
    o.sector = std::move(sector);
    o.momentum = std::move(momentum);
    return o;
}
ModelSpectrumOptions spectrum_options(std::map<std::string, double> sector, std::vector<double> momentum = {},
                                      bool all_momenta = false, std::optional<int> count = {}) {
    ModelSpectrumOptions o;
    o.sector = std::move(sector);
    o.momentum = std::move(momentum);
    o.all_momenta = all_momenta;
    o.count = count;
    return o;
}
ProductOptions product_sector(std::map<std::string, double> sector) {
    ProductOptions o;
    o.sector = std::move(sector);
    return o;
}
// A charge-only basis of product dimension prod(dims) with local tridiagonal
// Hermitian terms, for exact dispatch-threshold checks.
std::pair<mdl::CompiledOperator, ProjectedBasis> local_problem(std::vector<int> const& dims) {
    mdl::CompiledOperator op{dims, {}};
    std::vector<SpectrumLocalSpace> spaces;
    double weight = 1.;
    for (int site = 0; site < static_cast<int>(dims.size()); ++site, weight *= .1234) {
        int const dim = dims[static_cast<std::size_t>(site)];
        std::vector<C> matrix(static_cast<std::size_t>(dim * dim));
        for (int state = 0; state < dim; ++state) {
            matrix[static_cast<std::size_t>(state * dim + state)] = weight * state;
            if (state + 1 < dim)
                matrix[static_cast<std::size_t>(state * dim + state + 1)]
                    = matrix[static_cast<std::size_t>((state + 1) * dim + state)] = .3 * weight;
        }
        op.terms.push_back({1., {{site, std::move(matrix)}}, {}});
        spaces.push_back({std::vector<double>(static_cast<std::size_t>(dim), 1.), {}});
    }
    return {op, make_projected_basis(spaces, {})};
}
} // namespace

TEST(ProjectedProblem, NativeModelOwnsInputsAndResolvesTheCompleteSector) {
    auto model = mdl::heisenberg(alps::lattice::chain(6));
    std::vector<double> resolved;
    for (auto const& momentum : model.lat().momenta()) {
        auto problem = prepare_projected(model, projected({{"Sz", 0.}}, momentum));
        EXPECT_EQ(problem.charges(), std::vector<std::string>{"Sz"});
        EXPECT_EQ(problem.basis().momentum(), momentum);
        auto result = problem.eigenstates(options(100));
        EXPECT_EQ(result.spectrum.method, "dense");
        EXPECT_EQ(result.spectrum.momentum, momentum);
        EXPECT_EQ(result.eigenvectors_column_major.size(), problem.dimension() * result.spectrum.eigenvalues.size());
        resolved.insert(resolved.end(), result.spectrum.eigenvalues.begin(), result.spectrum.eigenvalues.end());
    }
    auto product = prepare_product(model, product_sector({{"Sz", 0.}})).eigenstates(options(100));
    ASSERT_EQ(resolved.size(), product.eigenvalues.size());
    resolved = sorted(resolved);
    for (std::size_t i = 0; i < resolved.size(); ++i) EXPECT_NEAR(resolved[i], product.eigenvalues[i], 1e-12);

    auto problem = prepare_projected(model, projected({{"Sz", 0.}}, {pi}));
    auto const before = problem.eigenstates(options(2)).spectrum.eigenvalues;
    model.add(100. * mdl::TermSum::located(mdl::spin(), "Sz", 0));
    EXPECT_EQ(problem.eigenstates(options(2)).spectrum.eigenvalues, before);
    auto copy = problem;
    EXPECT_EQ(copy.basis().representatives().data(), problem.basis().representatives().data());
    auto moved = std::move(problem);
    EXPECT_THROW(problem.eigenstates(), std::logic_error);
    EXPECT_EQ(moved.eigenstates(options(2)).spectrum.eigenvalues, before);
}

TEST(ProjectedProblem, CompiledModelChecksTheFullOperatorBeforeEmptySectorSelection) {
    ProjectedSpaces spaces{{"N"}, {{{1., -1.}, {{0., 1.}}}}};
    mdl::CompiledOperator op{{2}, {{1., {{0, {0., 1., 1., 0.}}}, {}}}};
    EXPECT_THROW(prepare_projected(op, spaces, {}, projected({{"N", 9.}})), std::invalid_argument);
    EXPECT_THROW(model_spectrum(op, spaces, {}, spectrum_options({{"N", 9.}})), std::invalid_argument);
    // The input parity is rightmost: X*F is anti-Hermitian although X is not.
    op.terms[0].jordan_wigner = {0};
    spaces.charges.clear();
    spaces.local_spaces[0].quantum_numbers.clear();
    EXPECT_THROW(prepare_projected(op, spaces, {}), std::invalid_argument);
    op.terms[0].coefficient = C(0., 1.);
    auto valid = prepare_projected(op, spaces, {}).eigenstates(options(2));
    EXPECT_NEAR(valid.spectrum.eigenvalues[0], -1., 1e-12);
    EXPECT_NEAR(valid.spectrum.eigenvalues[1], 1., 1e-12);
    op.terms[0].factors[0].site = 1;
    EXPECT_THROW(prepare_projected(op, spaces, {}), std::invalid_argument);
}

TEST(ProjectedProblem, CompiledAndExplicitModelSpectraHaveTheSameContract) {
    auto model = mdl::heisenberg(alps::lattice::chain(6));
    auto compiled = mdl::compile_hamiltonian(model).op;
    auto spaces = projected_spaces(mdl::physical_basis(model));
    std::vector<SpectrumTranslation> translations{{{1, 2, 3, 4, 5, 0}, 6}};
    auto explicit_result = model_spectrum(model, spectrum_options({}, {}, true));
    auto compiled_result = model_spectrum(compiled, spaces, translations, spectrum_options({}, {}, true));
    ASSERT_EQ(compiled_result.blocks.size(), explicit_result.blocks.size());
    EXPECT_EQ(compiled_result.charges, explicit_result.charges);
    EXPECT_EQ(compiled_result.directions, explicit_result.directions);
    for (std::size_t i = 0; i < explicit_result.blocks.size(); ++i) {
        auto const& a = compiled_result.blocks[i];
        auto const& b = explicit_result.blocks[i];
        EXPECT_EQ(a.momentum, b.momentum);
        EXPECT_EQ(a.quantum_numbers, b.quantum_numbers);
        EXPECT_EQ(a.block_dimension, b.block_dimension);
        EXPECT_EQ(a.eigenvalues, b.eigenvalues);
    }
}

TEST(ProjectedProblem, LatticeTranslationsMatchCellOrderedGenerators) {
    struct Case {
        alps::lattice::Lattice lattice;
        std::vector<std::vector<std::size_t>> generators;
        std::size_t dimension;
    };
    // Generator arrays are the Python lattice's, pinned by the eigenstate lock.
    std::vector<Case> cases;
    cases.push_back({alps::lattice::honeycomb(2), {{4, 5, 6, 7, 0, 1, 2, 3}, {2, 3, 0, 1, 6, 7, 4, 5}}, 22});
    cases.push_back({alps::lattice::ladder(3), {{2, 3, 4, 5, 0, 1}}, 8});
    for (auto const& item : cases) {
        ASSERT_EQ(item.lattice.translation_generators(), item.generators);
        auto model = mdl::heisenberg(item.lattice);
        std::vector<double> const zero(item.generators.size(), 0.);
        auto native = prepare_projected(model, projected({{"Sz", 0.}}, zero));
        std::vector<SpectrumTranslation> translations;
        for (std::size_t i = 0; i < item.generators.size(); ++i)
            translations.push_back({item.generators[i],
                static_cast<std::size_t>(item.lattice.extent()[static_cast<std::size_t>(item.lattice.periodic_directions()[i])])});
        std::vector<mdl::SiteBasisPtr> sites;
        for (int site = 0; site < model.num_sites(); ++site) sites.push_back(model.site_for(site));
        auto explicit_translations = prepare_projected(model.expand(), sites, model.conserved_qns(),
            model.fermion_convention(), translations, projected({{"Sz", 0.}}, zero));
        EXPECT_EQ(native.dimension(), item.dimension);
        EXPECT_EQ(native.basis().representatives(), explicit_translations.basis().representatives());
    }
    auto open = mdl::heisenberg(alps::lattice::chain(4, 1., alps::lattice::Boundary::open));
    EXPECT_EQ(message([&] { prepare_projected(open, projected({}, {0.})); }),
              "momentum resolution requires at least one periodic lattice direction");
    EXPECT_EQ(prepare_projected(open, projected({{"Sz", 0.}})).dimension(), 6);
}

TEST(ProjectedProblem, MomentumSelectionFollowsTheGridConvention) {
    std::vector<std::size_t> const ring{4}, grid{3, 4};
    struct Case {
        double wavevector;
        std::size_t index;
        double momentum;
    };
    for (auto const& item : {Case{pi / 2 + 5e-11, 1, pi / 2}, Case{-pi / 2, 3, 3 * pi / 2}, Case{2 * pi, 0, 0.},
                             Case{-1e-20, 0, 0.}, Case{9 * pi, 2, pi}}) {
        std::vector<double> const input{item.wavevector};
        auto selected = select_momentum(ring, input);
        EXPECT_EQ(selected.index, item.index);
        EXPECT_EQ(selected.momentum, std::vector<double>{item.momentum});
    }
    std::vector<double> const both{2 * pi / 3, pi / 2};
    auto selected = select_momentum(grid, both);
    EXPECT_EQ(selected.index, 5);
    EXPECT_EQ(selected.momentum, (std::vector<double>{2 * pi / 3, pi / 2}));
    EXPECT_EQ(select_momentum({}, {}).index, 0);
    std::vector<double> const off{pi / 2 + 2e-10}, nan{std::numeric_limits<double>::quiet_NaN()};
    EXPECT_EQ(message([&] { select_momentum(ring, off); }), "momentum is not on the lattice's allowed translation grid");
    EXPECT_EQ(message([&] { select_momentum(ring, nan); }), "momentum values must be finite real numbers");
    EXPECT_EQ(message([&] { select_momentum(grid, off); }), "momentum needs one component per periodic lattice direction");
    std::vector<std::size_t> const zero{0}, huge{std::size_t{1} << 53},
        overflowing{std::size_t{1} << 20, std::size_t{1} << 20, std::size_t{1} << 30};
    EXPECT_EQ(message([&] { select_momentum(zero, off); }), "translation periods must be positive");
    std::vector<double> const one{0.}, three{0., 0., 0.};
    EXPECT_EQ(message([&] { select_momentum(huge, one); }), "translation group is too large");
    EXPECT_EQ(message([&] { select_momentum(overflowing, three); }), "translation group is too large");
}

TEST(ProjectedProblem, ProjectedSpacesSortChargesAndFollowTheConvention) {
    auto site = mdl::spinful_fermion();
    for (auto convention : {mdl::FermionConvention::strict, mdl::FermionConvention::lattice}) {
        auto spaces = projected_spaces(mdl::physical_basis({site, site}, {"Sz", "N"}, convention));
        EXPECT_EQ(spaces.charges, (std::vector<std::string>{"N", "Sz"}));
        ASSERT_EQ(spaces.local_spaces.size(), 2);
        auto const& local = spaces.local_spaces[0];
        ASSERT_EQ(local.quantum_numbers.size(), 2);
        EXPECT_EQ(local.quantum_numbers[0], site->quantum_numbers.at("N"));
        EXPECT_EQ(local.quantum_numbers[1], site->quantum_numbers.at("Sz"));
        EXPECT_EQ(local.parity, convention == mdl::FermionConvention::strict
                                    ? site->parity_operator : std::vector<double>(static_cast<std::size_t>(site->dim), 1.));
    }
}

TEST(ProjectedProblem, DispatchSharesTheProductEigenstatePolicy) {
    for (auto const& item : {std::tuple<std::vector<int>, std::size_t, bool>{{10, 10, 10}, 1000, false},
                             std::tuple<std::vector<int>, std::size_t, bool>{{7, 11, 13}, 1001, true}}) {
        auto const& dims = std::get<0>(item);
        auto const dimension = std::get<1>(item);
        auto const lanczos = std::get<2>(item);
        auto [op, basis] = local_problem(dims);
        auto problem = prepare_projected(op, basis);
        ASSERT_EQ(problem.dimension(), dimension);
        EXPECT_TRUE(problem.charges().empty());
        auto result = problem.eigenstates(options(1));
        EXPECT_EQ(result.sparse_requested, lanczos);
        EXPECT_EQ(result.spectrum.method, lanczos ? "lanczos" : "dense");
        auto dense = problem.eigenstates(options(1, EigenstateMethod::dense));
        EXPECT_NEAR(result.spectrum.eigenvalues.at(0), dense.spectrum.eigenvalues.at(0), 1e-12);
        if (lanczos) {
            EXPECT_EQ(result.spectrum.tolerance, 1e-10);
            EXPECT_EQ(result.spectrum.max_iterations, 1000);
            EXPECT_EQ(result.residual_norms.size(), 1);
            auto explicit_controls = problem.eigenstates(options(1, EigenstateMethod::lanczos, "SA", 1e-9, 800));
            EXPECT_EQ(explicit_controls.spectrum.tolerance, 1e-9);
            EXPECT_EQ(explicit_controls.spectrum.max_iterations, 800);
        } else {
            EXPECT_FALSE(result.spectrum.tolerance);
            EXPECT_TRUE(result.residual_norms.empty());
        }
    }
    auto [op, basis] = local_problem({2, 2});
    auto tiny = prepare_projected(op, basis).eigenstates(options(3, EigenstateMethod::lanczos));
    EXPECT_TRUE(tiny.sparse_requested);
    EXPECT_EQ(tiny.spectrum.method, "dense");
    EXPECT_EQ(tiny.spectrum.nconv, 3);
    EXPECT_FALSE(tiny.spectrum.tolerance);
}

TEST(ProjectedProblem, RequestsPrecedeModelValidation) {
    auto spin = mdl::spin();
    std::vector<mdl::SiteBasisPtr> sites(4, spin);
    std::vector<SpectrumTranslation> ring{{{1, 2, 3, 0}, 4}};
    mdl::TermSum imaginary;
    for (int site = 0; site < 4; ++site) imaginary = imaginary + C(0., 1.) * mdl::TermSum::located(spin, "Sz", site);
    auto prepare = [&](mdl::TermSum const& terms, ProjectedOptions const& o) {
        return message([&] { prepare_projected(terms, sites, {"Sz"}, mdl::FermionConvention::strict, ring, o); });
    };
    EXPECT_EQ(prepare(imaginary, projected({{"N", 1.}}, {0.})), "sector contains a charge not declared conserved: N");
    EXPECT_EQ(prepare(imaginary, projected({}, {.1})), "momentum is not on the lattice's allowed translation grid");
    EXPECT_EQ(prepare(imaginary, projected({}, {0.})), "Hamiltonian must be Hermitian");
    mdl::TermSum flip;
    for (int site = 0; site < 4; ++site) flip = flip + mdl::TermSum::located(spin, "Sx", site);
    // Every declared charge is validated, even without a selected sector.
    EXPECT_EQ(prepare(flip, projected({}, {0.})), "Hamiltonian does not conserve declared quantum number Sz");
    EXPECT_NO_THROW(prepare_product(flip, sites, {"Sz"}));
    mdl::TermSum nonuniform = mdl::TermSum::located(spin, "Sz", 0);
    EXPECT_EQ(prepare(nonuniform, projected({}, {0.})),
              "projected ED: Hamiltonian is not invariant under the lattice translations");
}

TEST(ProjectedProblem, EmptyBlocksFailedSolvesAndConcurrentReuse) {
    auto model = mdl::heisenberg(alps::lattice::chain(4));
    auto empty = prepare_projected(model, projected({{"Sz", 2.}}, {pi / 2}));
    EXPECT_EQ(empty.dimension(), 0);
    EXPECT_EQ(empty.basis().parent_dimension(), 1);
    auto nothing = empty.eigenstates();
    EXPECT_TRUE(nothing.spectrum.eigenvalues.empty());
    EXPECT_EQ(nothing.spectrum.method, "dense");
    EXPECT_EQ(nothing.spectrum.nconv, 0);
    EXPECT_FALSE(nothing.sparse_requested);

    // The k = pi singlet block of the 16-site chain (about 800 states) is
    // above the Lanczos dispatch threshold; two Lanczos steps cannot converge.
    auto problem = prepare_projected(mdl::heisenberg(alps::lattice::chain(16)), projected({{"Sz", 0.}}, {pi}));
    ASSERT_GT(problem.dimension(), 8 * 3 + 200);
    EXPECT_THROW(problem.eigenstates(options(0)), std::invalid_argument);
    EXPECT_THROW(problem.eigenstates(options(6, EigenstateMethod::automatic, "invalid")), std::invalid_argument);
    EXPECT_THROW(problem.eigenstates(options(6, static_cast<EigenstateMethod>(-1))), std::invalid_argument);
    EXPECT_THROW(problem.eigenstates(options(3, EigenstateMethod::lanczos, "SA", {}, 1)), std::runtime_error);
    auto const reference = problem.eigenstates(options(3, EigenstateMethod::lanczos));
    EXPECT_EQ(reference.spectrum.block_dimension, problem.dimension());
    EXPECT_EQ(reference.spectrum.method, "lanczos");
    auto const dense = problem.eigenstates(options(3, EigenstateMethod::dense));
    for (std::size_t i = 0; i < 3; ++i) EXPECT_NEAR(reference.spectrum.eigenvalues[i], dense.spectrum.eigenvalues[i], 1e-8);
    std::vector<std::future<std::vector<double>>> solves;
    for (int i = 0; i < 4; ++i)
        solves.push_back(std::async(std::launch::async, [&, i] {
            auto method = i % 2 ? EigenstateMethod::dense : EigenstateMethod::lanczos;
            return problem.eigenstates(options(3, method)).spectrum.eigenvalues;
        }));
    for (int i = 0; i < 4; ++i)
        EXPECT_EQ(solves[static_cast<std::size_t>(i)].get(), (i % 2 ? dense : reference).spectrum.eigenvalues);
}

TEST(ModelSpectrum, NativeModelResolvesChargesAndMomenta) {
    auto model = mdl::heisenberg(alps::lattice::chain(6));
    auto product = prepare_product(model, product_sector({{"Sz", 0.}})).eigenstates(options(100));
    auto all = model_spectrum(model, spectrum_options({{"Sz", 0.}}, {}, true));
    EXPECT_EQ(all.charges, std::vector<std::string>{"Sz"});
    EXPECT_TRUE(all.momentum.empty());
    ASSERT_EQ(all.blocks.size(), 6);
    std::vector<double> resolved;
    for (auto const& block : all.blocks) {
        EXPECT_EQ(block.quantum_numbers, std::vector<double>{0.});
        resolved.insert(resolved.end(), block.eigenvalues.begin(), block.eigenvalues.end());
    }
    resolved = sorted(resolved);
    ASSERT_EQ(resolved.size(), product.eigenvalues.size());
    for (std::size_t i = 0; i < resolved.size(); ++i) EXPECT_NEAR(resolved[i], product.eigenvalues[i], 1e-12);

    // One momentum reports its canonical value and matches the prepared block.
    auto one = model_spectrum(model, spectrum_options({{"Sz", 0.}}, {-pi / 3}));
    ASSERT_EQ(one.blocks.size(), 1);
    EXPECT_EQ(one.momentum, std::vector<double>{5 * pi / 3});
    EXPECT_EQ(one.blocks[0].momentum, one.momentum);
    auto prepared = prepare_projected(model, projected({{"Sz", 0.}}, {-pi / 3}))
                        .eigenstates(options(100, EigenstateMethod::dense));
    ASSERT_EQ(one.blocks[0].eigenvalues.size(), prepared.spectrum.eigenvalues.size());
    for (std::size_t i = 0; i < prepared.spectrum.eigenvalues.size(); ++i)
        EXPECT_NEAR(one.blocks[0].eigenvalues[i], prepared.spectrum.eigenvalues[i], 1e-12);

    // Without momentum every declared charge is resolved, unrestricted by default.
    auto charges = model_spectrum(model, spectrum_options({}, {}, false, 1));
    ASSERT_EQ(charges.blocks.size(), 7);
    EXPECT_TRUE(charges.blocks[0].momentum.empty());
    EXPECT_NEAR(charges.blocks[3].eigenvalues.at(0), product.eigenvalues.at(0), 1e-12);
}

TEST(ModelSpectrum, RequestsPrecedeModelValidation) {
    auto spin = mdl::spin();
    std::vector<mdl::SiteBasisPtr> sites(4, spin);
    std::vector<SpectrumTranslation> ring{{{1, 2, 3, 0}, 4}};
    mdl::TermSum flip;
    for (int site = 0; site < 4; ++site) flip = flip + mdl::TermSum::located(spin, "Sx", site);
    auto spectrum = [&](std::vector<SpectrumTranslation> const& translations, ModelSpectrumOptions const& o) {
        return message([&] { model_spectrum(flip, sites, {"Sz"}, mdl::FermionConvention::strict, translations, o); });
    };
    EXPECT_EQ(spectrum(ring, spectrum_options({}, {}, true, 0)), "eigenvalue count must be positive");
    EXPECT_EQ(spectrum(ring, spectrum_options({{"N", 1.}}, {}, true)), "sector contains a charge not declared conserved: N");
    EXPECT_EQ(spectrum(ring, spectrum_options({}, {0.}, true)), "request either one momentum or all momenta");
    EXPECT_EQ(spectrum(ring, spectrum_options({})), "translations require a momentum or all-momenta request");
    EXPECT_EQ(spectrum(ring, spectrum_options({}, {.1})), "momentum is not on the lattice's allowed translation grid");
    EXPECT_EQ(spectrum({}, spectrum_options({}, {0.})), "momentum needs one component per periodic lattice direction");
    // Every declared charge is validated, even without a selected sector.
    EXPECT_EQ(spectrum({}, spectrum_options({})), "Hamiltonian does not conserve declared quantum number Sz");
    auto open = mdl::heisenberg(alps::lattice::chain(4, 1., alps::lattice::Boundary::open));
    EXPECT_EQ(message([&] { model_spectrum(open, spectrum_options({}, {}, true)); }),
              "momentum resolution requires at least one periodic lattice direction");
}
