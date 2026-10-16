// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// XML catalog models through ALPS::classic_dmrg against ALPS::ed on the same
// compiled operator and physical basis: the operators the dmrg application
// solves. The cases carry the signs a catalog operator can get wrong in the
// engine: wrap bonds whose factors arrive in descending site order, fermion
// rings with an even particle number, the catalog's two fermionic charges
// (Nup and Ndown, no N), ladders whose bonds skip sites of the engine's chain,
// projector-form bond terms of spin 1, and correlations with i > j.
#include <alps/dmrg/classic.hpp>
#include <alps/ed/projected_basis.hpp>
#include <alps/ed/projected_problem.hpp>
#include <alps/model/catalog.hpp>
#include <gtest/gtest.h>

#include <complex>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace cd = alps::dmrg::classic;
namespace mdl = alps::models;
using Catalog = mdl::CatalogModel<>;
using Sector = std::map<std::string, double>;
using Entries = std::vector<std::pair<std::string, std::string>>;

alps::Parameters parameters(Entries const& entries) {
    alps::Parameters p;
    for (auto const& [key, value] : entries) p[key] = value;
    return p;
}

struct Reference {
    alps::ed::ProjectedEigenResult states;
    std::size_t dimension;
    double value(mdl::CompiledOperator const& op, std::size_t state) const {
        std::vector<std::complex<double>> column(states.eigenvectors_column_major.begin() + state * dimension,
                                                 states.eigenvectors_column_major.begin() + (state + 1) * dimension);
        return alps::ed::expect_projected_operator(op, states.problem.basis(), column).real();
    }
};

Reference exact(mdl::CompiledOperator const& op, mdl::PhysicalBasis const& basis, Sector const& sector, int count) {
    alps::ed::ProjectedOptions options;
    options.sector = sector;
    auto problem = alps::ed::prepare_projected(op, alps::ed::projected_spaces(basis), {}, options);
    alps::ed::EigenstateOptions request;
    request.count = count;
    request.method = alps::ed::EigenstateMethod::dense;
    auto const dimension = problem.dimension();
    return Reference{problem.eigenstates(request), dimension};
}

cd::Options sweeps(std::size_t states, Sector sector, std::size_t count = 1) {
    cd::Options options;
    options.sweeps = cd::fixed_sweeps(4, states);
    options.warmup_states = states;
    options.num_eigenvalues = count;
    options.sector = std::move(sector);
    options.lanczos_tolerance = 1e-13;
    return options;
}

// Ground state (and further targets) of a catalog model by DMRG and by ED.
void expect_energies(Entries const& entries, std::vector<std::string> const& conserved, Sector const& sector,
                     std::size_t states, std::size_t count = 1) {
    Catalog model(parameters(entries));
    auto const op = model.compile_hamiltonian();
    auto const basis = model.physical_basis(conserved);
    auto const reference = exact(op, basis, sector, static_cast<int>(count));
    auto const result = cd::solve(op, basis, sweeps(states, sector, count));
    ASSERT_EQ(result.energies.size(), count);
    for (std::size_t k = 0; k < count; ++k)
        EXPECT_NEAR(result.energies[k], reference.states.spectrum.eigenvalues[k], 1e-8) << "target " << k;
}

Entries spin_chain(std::string lattice, int length, std::string spin = "0.5") {
    return {{"MODEL", "spin"}, {"LATTICE", std::move(lattice)}, {"L", std::to_string(length)}, {"J", "1"},
            {"local_S", std::move(spin)}};
}

TEST(CatalogClassicDmrg, PeriodicSpinChainWrapBondArrivesDescending) {
    expect_energies(spin_chain("chain lattice", 8), {"Sz"}, {{"Sz", 0.}}, 64);
}

TEST(CatalogClassicDmrg, SpinOneProjectorTermsCombine) {
    expect_energies(spin_chain("open chain lattice", 8, "1"), {"Sz"}, {{"Sz", 0.}}, 200, 2);
}

TEST(CatalogClassicDmrg, OpenSpinLadderBondsSkipSites) {
    expect_energies({{"MODEL", "spin"}, {"LATTICE", "open ladder"}, {"L", "4"}, {"J0", "1"}, {"J1", "0.7"}},
                    {"Sz"}, {{"Sz", 0.}}, 128);
}

TEST(CatalogClassicDmrg, SpinlessFermionRingWithEvenParticleNumber) {
    Entries const ring{{"MODEL", "spinless fermions"}, {"LATTICE", "chain lattice"}, {"L", "8"}, {"t", "1"},
                       {"V", "0.5"}};
    expect_energies(ring, {"N"}, {{"N", 4.}}, 64);
    expect_energies(ring, {"N"}, {{"N", 3.}}, 64);
}

TEST(CatalogClassicDmrg, HubbardCarriesTwoFermionicCharges) {
    Entries const open{{"MODEL", "fermion Hubbard"}, {"LATTICE", "open chain lattice"}, {"L", "6"}, {"t", "1"},
                       {"U", "4"}};
    expect_energies(open, {"Ndown", "Nup"}, {{"Nup", 3.}, {"Ndown", 3.}}, 200, 2);
    expect_energies(open, {"Ndown", "Nup"}, {{"Nup", 3.}, {"Ndown", 2.}}, 200);
    Entries ring = open;
    ring[1].second = "chain lattice";
    expect_energies(ring, {"Ndown", "Nup"}, {{"Nup", 2.}, {"Ndown", 2.}}, 200);
}

TEST(CatalogClassicDmrg, HubbardLadderStringsSpanSkippedSites) {
    expect_energies({{"MODEL", "fermion Hubbard"}, {"LATTICE", "open ladder"}, {"L", "3"}, {"t0", "1"},
                     {"t1", "0.8"}, {"U", "3"}},
                    {"Ndown", "Nup"}, {{"Nup", 3.}, {"Ndown", 3.}}, 256);
}

TEST(CatalogClassicDmrg, TJChainConservesNAndSz) {
    expect_energies({{"MODEL", "t-J"}, {"LATTICE", "open chain lattice"}, {"L", "8"}, {"t", "1"}, {"J", "0.5"}},
                    {"N", "Sz"}, {{"N", 6.}, {"Sz", 0.}}, 200);
}

TEST(CatalogClassicDmrg, CorrelationsInBothSiteOrders) {
    // One expectation value per ordered pair, compiled as the application
    // compiles them: a two-site bond operator on (i, j), descending for i > j.
    Catalog model(parameters({{"MODEL", "spinless fermions"}, {"LATTICE", "open chain lattice"}, {"L", "6"},
                              {"t", "1"}, {"V", "0.5"}}));
    auto const op = model.compile_hamiltonian();
    auto const basis = model.physical_basis({"N"});
    Sector const sector{{"N", 3.}};
    auto const reference = exact(op, basis, sector, 1);
    cd::Measurements requests;
    std::vector<std::pair<int, int>> const pairs{{1, 4}, {4, 1}, {0, 5}, {5, 0}, {2, 2}};
    for (auto const& [i, j] : pairs) {
        auto const name = std::to_string(i) + "," + std::to_string(j);
        if (i == j) {
            requests.expectation["n " + name] =
                model.compile_operator(alps::SiteOperator("cdag(i)*c(i)", "i"), Catalog::site_descriptor(i));
        } else {
            requests.expectation["hop " + name] = model.compile_operator(
                alps::BondOperator("cdag(i)*c(j)", "i", "j"), Catalog::site_descriptor(i), Catalog::site_descriptor(j));
            requests.expectation["nn " + name] = model.compile_operator(
                alps::BondOperator("n(i)*n(j)", "i", "j"), Catalog::site_descriptor(i), Catalog::site_descriptor(j));
        }
    }
    auto const result = cd::solve(op, basis, sweeps(64, sector), requests);
    EXPECT_NEAR(result.energies.at(0), reference.states.spectrum.eigenvalues.at(0), 1e-8);
    for (auto const& [name, request] : requests.expectation) {
        auto const& measured = result.measurements.at(name);
        ASSERT_EQ(measured.values.size(), 1u);
        EXPECT_NEAR(measured.values[0].at(0), reference.value(request, 0), 1e-7) << name;
    }
    // <c+_1 c_4> is real here, so the two orders agree.
    EXPECT_NEAR(result.measurements.at("hop 1,4").values[0][0], result.measurements.at("hop 4,1").values[0][0], 1e-7);
}

TEST(CatalogClassicDmrg, CatalogPhysicalBasisFollowsTheSiteBasis) {
    Catalog model(parameters({{"MODEL", "fermion Hubbard"}, {"LATTICE", "open chain lattice"}, {"L", "4"}}));
    auto const basis = model.physical_basis({"Nup", "Ndown", "Nup"});
    EXPECT_EQ(basis.conserved_quantum_numbers, (std::vector<std::string>{"Ndown", "Nup"}));
    ASSERT_EQ(basis.sites.size(), 4u);
    auto const& site = basis.sites[0];
    EXPECT_EQ(site.statistics, mdl::Statistics::fermionic);
    EXPECT_EQ(site.parity, (std::vector<double>{1., -1., -1., 1.}));
    EXPECT_EQ(site.quantum_numbers.at("Nup").size(), 4u);
    EXPECT_EQ(site.labels.front(), "Nup=0 Ndown=0");
    // A charge the site basis lacks is carried as zeros.
    auto const spin = Catalog(parameters(spin_chain("open chain lattice", 4))).physical_basis({"N", "Sz"});
    EXPECT_EQ(spin.sites[0].quantum_numbers.at("N"), (std::vector<double>{0., 0.}));
    EXPECT_EQ(spin.sites[0].statistics, mdl::Statistics::bosonic);
}

} // namespace
