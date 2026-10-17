// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/model/tests/models/catalog_model_parity.cpp and its frozen
// params_*.fixture inputs (captured from Python before the model-export
// metadata was removed). Adaptation: GoogleTest; the catalog ground state is
// obtained through the catalog importer and ALPS::ed instead of alpsx's
// Eigen-based CatalogMatrix; the lattice half of each catalog input comes
// from the kit lattice's catalog name and extent (alpsx used the removed
// to_alps_params export); the shipped models.xml/lattices.xml replace the
// TOML catalogs.
//
// Numerical lock for the C++ operator models against independent, frozen
// catalog descriptions. Covers boundary conditions, parallel periodic bonds,
// restricted bond types and distinct leg/rung interactions. The catalog and
// explicit local-operator assembly must produce the same ground-state energy.
#include <alps/ed/eigensolver.hpp>
#include <alps/ed/product_basis.hpp>
#include <alps/ed/product_problem.hpp>
#include <alps/ed/projected_problem.hpp>
#include <alps/ed/sparse_matrix.hpp>
#include <alps/model/catalog.hpp>
#include <alps/models/builtin.hpp>
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace lat = alps::lattice;
namespace mdl = alps::models;

namespace {
using Params = std::vector<std::pair<std::string, std::string>>;
struct Case {
    std::string name;
    std::function<mdl::Hamiltonian()> factory;
    Params catalog; // the frozen params_<name>.fixture, verbatim
};

std::vector<Case> const& cases() {
    using B = lat::Boundary;
    static std::vector<Case> const r = {
        {"heisenberg_chain_L6", [] { return mdl::heisenberg(lat::chain(6), 1.0); },
         {{"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"J", "1"}, {"MODEL", "spin"}, {"local_S", "0.5"}}},
        {"heisenberg_ladder_L3_p", [] { return mdl::heisenberg(lat::ladder(3), 1.0); },
         {{"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"J0", "1"}, {"J1", "1"}, {"MODEL", "spin"}, {"local_S", "0.5"}}},
        {"heisenberg_ladder_L3_o", [] { return mdl::heisenberg(lat::ladder(3, 1.0, B::open), 1.0); },
         {{"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"J0", "1"}, {"J1", "1"}, {"MODEL", "spin"}, {"local_S", "0.5"}}},
        // A periodic extent-2 direction contains both wrap images. Dropping
        // either parallel bond halves that direction's interaction.
        {"heisenberg_square_2x3_p", [] { return mdl::heisenberg(lat::square(2, 3), 1.0); },
         {{"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"J", "1"}, {"MODEL", "spin"}, {"local_S", "0.5"}}},
        {"heisenberg_square_3x3_p", [] { return mdl::heisenberg(lat::square(3, 3), 1.0); },
         {{"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"J", "1"}, {"MODEL", "spin"}, {"local_S", "0.5"}}},
        // The next-nearest bonds of a periodic four-site chain also have
        // two images for each site pair.
        {"j1j2_nnn_L4", [] { return mdl::j1j2(lat::nnn_chain(4), 1.0, 0.5); },
         {{"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"J", "1"}, {"J'", "0.5"}, {"MODEL", "spin"}, {"local_S", "0.5"}}},
        {"heisenberg_square_2x3_o", [] { return mdl::heisenberg(lat::square(2, 3, 1.0, B::open), 1.0); },
         {{"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"J", "1"}, {"MODEL", "spin"}, {"local_S", "0.5"}}},
        {"heisenberg_ladder_rungs_only", [] { return mdl::heisenberg(lat::ladder(3), 0.7, 0.0, 0.5, 1); },
         {{"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"J1", "0.69999999999999996"}, {"MODEL", "spin"}, {"local_S", "0.5"}}},
        {"ising_ladder_L3_h", [] { return mdl::ising(lat::ladder(3), 1.0, 0.3); },
         {{"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"Jxy0", "0"}, {"Jxy1", "0"}, {"Jz0", "1"}, {"Jz1", "1"},
          {"MODEL", "spin"}, {"h", "0.29999999999999999"}, {"local_S", "0.5"}}},
        {"xxz_ladder_L3", [] { return mdl::xxz(lat::ladder(3), 1.0, 0.5); },
         {{"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"Jxy0", "1"}, {"Jxy1", "1"}, {"Jz0", "0.5"}, {"Jz1", "0.5"},
          {"MODEL", "spin"}, {"local_S", "0.5"}}},
        {"tfi_ladder_L3", [] { return mdl::transverse_field_ising(lat::ladder(3), 1.0, 0.7); },
         {{"Gamma", "0.69999999999999996"}, {"Jxy0", "0"}, {"Jxy1", "0"}, {"Jz0", "1"}, {"Jz1", "1"},
          {"MODEL", "spin"}, {"local_S", "0.5"}}},
        {"heisenberg_d_chain_L4_S1", [] { return mdl::heisenberg_d(lat::chain(4), 1.0, 0.2, 0.0, 1.0); },
         {{"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"D", "0.20000000000000001"}, {"J", "1"}, {"MODEL", "spin"},
          {"local_S", "1"}}},
        {"j1j2_nnn_L6", [] { return mdl::j1j2(lat::nnn_chain(6), 1.0, 0.5); },
         {{"CONSERVED_QUANTUMNUMBERS", "Sz"}, {"J", "1"}, {"J'", "0.5"}, {"MODEL", "spin"}, {"local_S", "0.5"}}},
    };
    return r;
}

// The lattice half of a catalog input: the kit lattice's catalog name and
// extent. alps_name() names the periodic graph for every factory but the
// chains, so open ladders and squares select the catalog's open variant
// (the rung direction of a ladder is always open; its first direction decides).
alps::Parameters lattice_parameters(lat::Lattice const& lattice) {
    alps::Parameters p;
    std::string name = lattice.alps_name();
    if (lattice.boundary().at(0) == lat::Boundary::open
        && (name == "ladder" || name == "square lattice" || name == "simple cubic lattice"))
        name = "open " + name;
    p["LATTICE"] = name;
    p["L"] = lattice.extent().at(0);
    if (lattice.extent().size() > 1) p["W"] = lattice.extent().at(1);
    if (lattice.extent().size() > 2) p["H"] = lattice.extent().at(2);
    return p;
}

double kit_e0(mdl::Hamiltonian const& model) {
    alps::ed::EigenstateOptions options;
    options.count = 1;
    return alps::ed::prepare_product(model).eigenstates(options).eigenvalues.at(0);
}

// Lowest eigenvalue of the catalog operator on the full product space.
double catalog_e0(alps::Parameters const& p) {
    mdl::CatalogModel<> model(p);
    auto const op = model.compile_hamiltonian();
    auto const basis = alps::ed::enumerate_product_basis(op.site_dims);
    std::vector<std::int64_t> const states(basis.states.begin(), basis.states.end());
    auto const matrix = alps::ed::build_csr_from_local<std::complex<double>>(
        op, states, basis.dimension, model.local_parity());
    std::vector<std::complex<double>> dense(basis.dimension * basis.dimension);
    for (std::size_t row = 0; row < matrix.rows; ++row)
        for (auto i = matrix.offsets[row]; i < matrix.offsets[row + 1]; ++i)
            dense[static_cast<std::size_t>(matrix.indices[static_cast<std::size_t>(i)]) * basis.dimension + row] =
                matrix.values[static_cast<std::size_t>(i)];
    return alps::ed::dense_eigensystem(std::move(dense), basis.dimension, false, 1, "SA").eigenvalues.at(0);
}
} // namespace

TEST(CatalogModelParity, FrozenCatalogInputsReproduceTheKitGroundStates) {
    for (auto const& item : cases()) {
        SCOPED_TRACE(item.name);
        mdl::Hamiltonian const model = item.factory();
        auto p = lattice_parameters(model.lat());
        for (auto const& [key, value] : item.catalog) p[key] = value;
        // Regression hints: catalog above kit = dropped couplings, below =
        // phantom periodic bonds.
        EXPECT_NEAR(catalog_e0(p), kit_e0(model), 1e-9);
    }
}

TEST(CatalogModelParity, TypedCouplingsSelectDistinctLegAndRungInteractions) {
    // Distinct leg/rung XXZ interactions against explicit catalog inputs.
    lat::Lattice const ladder = lat::ladder(3);
    mdl::SiteBasisPtr const site = mdl::spin(0.5);
    auto xxz_bond = [site](double jz, double jxy) {
        return [site, jz, jxy](int i, int j) {
            using mdl::TermSum;
            return TermSum::located(site, "Sz", i) * TermSum::located(site, "Sz", j) * std::complex<double>(jz)
                + TermSum::located(site, "Sp", i) * TermSum::located(site, "Sm", j) * std::complex<double>(0.5 * jxy)
                + TermSum::located(site, "Sm", i) * TermSum::located(site, "Sp", j) * std::complex<double>(0.5 * jxy);
        };
    };
    mdl::Hamiltonian model(ladder, site, mdl::FermionConvention::strict, {"Sz"});
    model.add_bond(1.0, xxz_bond(1.0, 0.6), 0); // legs
    model.add_bond(1.0, xxz_bond(0.4, 1.1), 1); // rungs
    auto p = lattice_parameters(ladder);
    p["MODEL"] = "spin";
    p["local_S"] = 0.5;
    p["Jz0"] = 1.0;
    p["Jxy0"] = 0.6;
    p["Jz1"] = 0.4;
    p["Jxy1"] = 1.1;
    p["h"] = 0.0;
    p["Gamma"] = 0.0;
    p["D"] = 0.0;
    EXPECT_NEAR(catalog_e0(p), kit_e0(model), 1e-9);
}

namespace {
std::vector<double> all_levels(alps::ed::ProductProblem const& problem) {
    alps::ed::EigenstateOptions options;
    options.count = static_cast<int>(problem.dimension());
    options.method = alps::ed::EigenstateMethod::dense;
    return problem.eigenstates(options).eigenvalues;
}
std::vector<double> all_levels(alps::ed::ProjectedProblem const& problem) {
    alps::ed::EigenstateOptions options;
    options.count = static_cast<int>(problem.dimension());
    options.method = alps::ed::EigenstateMethod::dense;
    return problem.eigenstates(options).spectrum.eigenvalues;
}
void expect_same_levels(std::vector<double> const& a, std::vector<double> const& b) {
    ASSERT_EQ(a.size(), b.size());
    for (std::size_t i = 0; i < a.size(); ++i) EXPECT_NEAR(a[i], b[i], 1e-10) << "level " << i;
}
} // namespace

// A compiled catalog operator (Jordan-Wigner strings included) in a charge
// sector: the product path, the projected path without translations and, for
// the Hubbard chain, the kit model in the same sector give one spectrum.
TEST(CatalogModelParity, CompiledCatalogSectorsAgreeAcrossProductProjectedAndKit) {
    struct Item { Params params; std::map<std::string, double> sector; };
    std::vector<Item> const items = {
        {{{"MODEL", "fermion Hubbard"}, {"LATTICE", "chain lattice"}, {"L", "4"}, {"t", "1"}, {"U", "4"},
          {"CONSERVED_QUANTUMNUMBERS", "Nup,Ndown"}, {"Nup_total", "2"}, {"Ndown_total", "2"}}, {{"Ndown", 2.}, {"Nup", 2.}}},
        {{{"MODEL", "t-J"}, {"LATTICE", "open chain lattice"}, {"L", "5"}, {"t", "1"}, {"J", "0.4"},
          {"CONSERVED_QUANTUMNUMBERS", "N,Sz"}, {"N_total", "3"}, {"Sz_total", "0.5"}}, {{"N", 3.}, {"Sz", 0.5}}},
    };
    for (auto const& item : items) {
        alps::Parameters p;
        for (auto const& [key, value] : item.params) p[key] = value;
        SCOPED_TRACE(static_cast<std::string>(p["MODEL"]));
        mdl::CatalogModel<> model(p);
        auto const op = model.compile_hamiltonian();
        auto const physical = model.physical_basis();
        auto const sector = model.sector();
        ASSERT_EQ(sector, item.sector);
        alps::ed::ProductOptions product;
        product.sector = sector;
        auto const levels = all_levels(alps::ed::prepare_product(op, physical, product));
        alps::ed::ProjectedOptions projected;
        projected.sector = sector;
        expect_same_levels(levels, all_levels(alps::ed::prepare_projected(op, alps::ed::projected_spaces(physical), {}, projected)));
        // The full space with no selected charge needs Hermiticity only.
        auto const full = alps::ed::prepare_product(op, physical);
        EXPECT_EQ(full.dimension(), static_cast<std::size_t>(std::pow(op.site_dims.front(), op.site_dims.size())));
    }
    // The kit Hubbard chain in N = 4, Sz = 0.
    mdl::CatalogModel<> catalog([] {
        alps::Parameters p;
        p["MODEL"] = "fermion Hubbard"; p["LATTICE"] = "chain lattice"; p["L"] = 4; p["t"] = 1; p["U"] = 4;
        p["CONSERVED_QUANTUMNUMBERS"] = "Nup,Ndown"; p["Nup_total"] = 2; p["Ndown_total"] = 2;
        return p;
    }());
    alps::ed::ProductOptions catalog_sector, kit_sector;
    catalog_sector.sector = catalog.sector();
    kit_sector.sector = {{"N", 4.}, {"Sz", 0.}};
    expect_same_levels(all_levels(alps::ed::prepare_product(catalog.compile_hamiltonian(), catalog.physical_basis(), catalog_sector)),
                       all_levels(alps::ed::prepare_product(mdl::hubbard(lat::chain(4), 1., 4.), kit_sector)));
    // Undeclared charges and mismatched site counts are rejected.
    alps::ed::ProductOptions wrong;
    wrong.sector = {{"Sz", 0.}};
    EXPECT_THROW(alps::ed::prepare_product(catalog.compile_hamiltonian(), catalog.physical_basis(), wrong), std::invalid_argument);
    auto shorter = catalog.compile_hamiltonian();
    shorter.site_dims.pop_back();
    EXPECT_THROW(alps::ed::prepare_product(shorter, catalog.physical_basis()), std::invalid_argument);
}
