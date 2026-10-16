// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The kit-lattice bridge (<alps/lattice/kit.hpp>) through the legacy lattice
// reader: every factory and two custom cells expand from their emitted
// <LATTICEGRAPH> document to the kit's sites, types, coordinates and ordered
// bonds; resolution selects a catalog graph only where the catalog's own
// expansion is identical; a written library drives graph_helper from
// parameters; the catalog importer runs on a kit lattice and reproduces the
// kit's ground state. The custom cells are those of alpsx's lattice fixtures
// (custom_bilayer_3x3_p, custom_cylinder_4x3).
#include <alps/ed/eigensolver.hpp>
#include <alps/ed/product_basis.hpp>
#include <alps/ed/product_problem.hpp>
#include <alps/ed/sparse_matrix.hpp>
#include <alps/lattice/kit.hpp>
#include <alps/model/catalog.hpp>
#include <alps/models/builtin.hpp>
#include <alps/models/lattice.hpp>
#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace lat = alps::lattice;
namespace mdl = alps::models;
using lat::Boundary;

namespace {
struct Case {
    std::string name;
    std::function<lat::Lattice()> factory;
    std::string catalog_graph; // "" when resolution must take the custom route
};

lat::Lattice custom_bilayer(int n) {
    lat::UnitCell uc(2);
    uc.add_vertex({0.0, 0.0}, 0);
    uc.add_vertex({0.0, 0.0}, 1);
    uc.add_edge(0, 0, 0, {}, {1, 0});
    uc.add_edge(0, 0, 0, {}, {0, 1});
    uc.add_edge(1, 1, 1, {}, {1, 0});
    uc.add_edge(1, 1, 1, {}, {0, 1});
    uc.add_edge(0, 1, 2);
    return uc.tile({n, n});
}

lat::Lattice custom_cylinder(int l, int w) {
    lat::UnitCell uc(2);
    uc.add_vertex({0.0, 0.0});
    uc.add_edge(0, 0, 0, {}, {1, 0});
    uc.add_edge(0, 0, 0, {}, {0, 1});
    return uc.tile({l, w}, {}, {Boundary::periodic, Boundary::open});
}

std::vector<Case> const& cases() {
    constexpr auto O = Boundary::open;
    static std::vector<Case> const r = {
        {"chain_L6_p", [] { return lat::chain(6); }, "chain lattice"},
        {"chain_L6_o", [] { return lat::chain(6, 1.0, O); }, "open chain lattice"},
        {"chain_L2_p", [] { return lat::chain(2); }, "chain lattice"},
        {"nnn_chain_L6_p", [] { return lat::nnn_chain(6); }, "nnn chain lattice"},
        {"nnn_chain_L5_o", [] { return lat::nnn_chain(5, 1.0, O); }, "nnn open chain lattice"},
        {"ladder_L4_p", [] { return lat::ladder(4); }, "ladder"},
        {"ladder_L4_o", [] { return lat::ladder(4, 1.0, O); }, "open ladder"},
        {"square_3x4_p", [] { return lat::square(3, 4); }, "square lattice"},
        {"square_3x4_o", [] { return lat::square(3, 4, 1.0, O); }, "open square lattice"},
        {"frustrated_square_3x3_p", [] { return lat::frustrated_square(3, 3); }, "frustrated square lattice"},
        {"frustrated_square_3x3_o", [] { return lat::frustrated_square(3, 3, 1.0, O); }, ""},
        {"triangular_3x3_p", [] { return lat::triangular(3, 3); }, "triangular lattice"},
        {"triangular_3x3_o", [] { return lat::triangular(3, 3, 1.0, O); }, ""},
        // The catalog leaves honeycomb and kagome vertices untyped and spells
        // "Kagome lattice": the kit's typed cells take the custom route.
        {"honeycomb_2x3_p", [] { return lat::honeycomb(2, 3); }, ""},
        {"honeycomb_2x2_o", [] { return lat::honeycomb(2, 2, 1.0, O); }, ""},
        {"kagome_2x2_p", [] { return lat::kagome(2, 2); }, ""},
        {"simple_cubic_2x3x2_p", [] { return lat::simple_cubic(2, 3, 2); }, "simple cubic lattice"},
        {"simple_cubic_2x2x2_o", [] { return lat::simple_cubic(2, 2, 2, 1.0, O); }, "open simple cubic lattice"},
        // No catalog graph of these names.
        {"rectangular_3x2_b05_p", [] { return lat::rectangular(3, 2, 1.0, 0.5); }, ""},
        {"fcc_2x2x2_p", [] { return lat::face_centered_cubic(2); }, ""},
        {"bcc_2x2x2_p", [] { return lat::body_centered_cubic(2); }, ""},
        {"tetragonal_2x2x2_c15_p", [] { return lat::tetragonal(2, 2, 2, 1.0, 1.5); }, ""},
        {"orthorhombic_2x2x2_abc_p", [] { return lat::orthorhombic(2, 2, 2, 1.0, 2.0, 3.0); }, ""},
        {"hexagonal_2x2x2_p", [] { return lat::hexagonal(2); }, ""},
        {"custom_bilayer_3x3_p", [] { return custom_bilayer(3); }, ""},
        {"custom_cylinder_4x3", [] { return custom_cylinder(4, 3); }, ""},
    };
    return r;
}

double kit_e0(mdl::Hamiltonian const& model) {
    alps::ed::EigenstateOptions options;
    options.count = 1;
    return alps::ed::prepare_product(model).eigenstates(options).eigenvalues.at(0);
}

// Lowest eigenvalue of a catalog operator on the full product space.
double catalog_e0(mdl::CatalogModel<> const& model) {
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

TEST(LatticeXmlBridge, EmittedDocumentsExpandToTheKitStructure) {
    for (auto const& item : cases()) {
        SCOPED_TRACE(item.name);
        auto const lattice = item.factory();
        auto const resolution = lat::custom_lattice(lattice);
        EXPECT_FALSE(resolution.from_catalog());
        EXPECT_EQ(resolution.graph, lat::custom_lattice_name(lattice));
        EXPECT_EQ(static_cast<std::string>(resolution.parameters["LATTICE"]), resolution.graph);
        auto const graph = lat::make_graph_helper(resolution);
        EXPECT_EQ(lat::structure_difference(*graph, lattice), "");
    }
}

TEST(LatticeXmlBridge, ResolutionUsesTheCatalogOnlyWhereItIsIdentical) {
    for (auto const& item : cases()) {
        SCOPED_TRACE(item.name);
        auto const lattice = item.factory();
        auto const resolution = lat::resolve_lattice(lattice);
        std::string catalog_difference = "(no catalog graph expected)";
        if (!item.catalog_graph.empty()) {
            alps::Parameters p;
            p["LATTICE"] = item.catalog_graph;
            p["L"] = lattice.extent().at(0);
            if (lattice.dimension() > 1) p["W"] = lattice.extent().at(1);
            if (lattice.dimension() > 2) p["H"] = lattice.extent().at(2);
            catalog_difference = lat::structure_difference(alps::graph_helper<>(p), lattice);
        }
        EXPECT_EQ(resolution.from_catalog(), !item.catalog_graph.empty()) << "catalog expansion: " << catalog_difference;
        if (resolution.from_catalog()) {
            EXPECT_EQ(resolution.graph, item.catalog_graph);
            EXPECT_EQ(static_cast<std::string>(resolution.parameters["LATTICE"]), item.catalog_graph);
            EXPECT_EQ(static_cast<int>(resolution.parameters["L"]), lattice.extent().at(0));
            if (lattice.dimension() > 1) EXPECT_EQ(static_cast<int>(resolution.parameters["W"]), lattice.extent().at(1));
            if (lattice.dimension() > 2) EXPECT_EQ(static_cast<int>(resolution.parameters["H"]), lattice.extent().at(2));
            for (int d = 0; d < lattice.dimension(); ++d)
                EXPECT_EQ(static_cast<std::string>(resolution.parameters["bc" + std::to_string(d + 1)]),
                          lat::to_string(lattice.boundary().at(static_cast<std::size_t>(d))));
            EXPECT_FALSE(resolution.parameters.defined("LATTICE_LIBRARY"));
        } else {
            EXPECT_EQ(resolution.graph, lat::custom_lattice_name(lattice));
            EXPECT_FALSE(resolution.parameters.defined("L"));
        }
        auto const graph = lat::make_graph_helper(resolution);
        EXPECT_EQ(lat::structure_difference(*graph, lattice), "");
    }
}

TEST(LatticeXmlBridge, MissingAndDivergentCatalogNamesTakeTheCustomRoute) {
    lat::UnitCell cell(1);
    cell.add_vertex({0.});
    cell.add_edge(0, 0, 0, {}, {1});
    auto const unknown = lat::Lattice::build(cell, {{1.}}, {4}, {Boundary::periodic}, "no such lattice");
    EXPECT_EQ(lat::catalog_graph_name(unknown), "no such lattice");
    auto const resolution = lat::resolve_lattice(unknown);
    EXPECT_FALSE(resolution.from_catalog());
    EXPECT_EQ(lat::structure_difference(*lat::make_graph_helper(resolution), unknown), "");
    // A periodic extent of 1: the kit drops the self-bond and so does the
    // legacy expansion of this tree, on both routes.
    auto const single = lat::chain(1);
    EXPECT_EQ(single.num_bonds(), 0);
    auto const fallback = lat::resolve_lattice(single);
    EXPECT_TRUE(fallback.from_catalog());
    EXPECT_EQ(lat::structure_difference(*lat::make_graph_helper(fallback), single), "");
    EXPECT_EQ(lat::structure_difference(*lat::make_graph_helper(lat::custom_lattice(single)), single), "");
    // A catalog whose library file is missing reports it rather than resolving.
    alps::Parameters missing;
    missing["LATTICE_LIBRARY"] = "/nonexistent/lattices.xml";
    EXPECT_THROW(lat::resolve_lattice(lat::chain(4), missing), std::runtime_error);
}

TEST(LatticeXmlBridge, WrittenLibraryDrivesTheLegacyReaderFromParameters) {
    auto const lattice = lat::honeycomb(2, 2);
    auto const resolution = lat::resolve_lattice(lattice);
    ASSERT_FALSE(resolution.from_catalog());
    auto const file = std::filesystem::path(ALPS_TEST_BINARY_DIR) / "kit_bridge_honeycomb_2x2.xml";
    std::filesystem::remove(file);
    auto const parameters = lat::write_lattice_library(resolution, file);
    ASSERT_TRUE(std::filesystem::exists(file));
    EXPECT_EQ(static_cast<std::string>(parameters["LATTICE_LIBRARY"]), std::filesystem::absolute(file).string());
    EXPECT_EQ(static_cast<std::string>(parameters["LATTICE"]), resolution.graph);
    alps::graph_helper<> const graph(parameters);
    EXPECT_EQ(lat::structure_difference(graph, lattice), "");
    // Vertex coordinates survive the file: the second sublattice is offset.
    EXPECT_NEAR(graph.coordinate(1)[0], lattice.sites()[1].coordinate[0], 1e-12);
    EXPECT_NEAR(graph.coordinate(1)[1], 1. / std::sqrt(3.), 1e-12);
    // The catalog route writes nothing.
    auto const named = lat::resolve_lattice(lat::square(3, 3));
    ASSERT_TRUE(named.from_catalog());
    auto const untouched = std::filesystem::path(ALPS_TEST_BINARY_DIR) / "kit_bridge_untouched.xml";
    std::filesystem::remove(untouched);
    auto const same = lat::write_lattice_library(named, untouched);
    EXPECT_FALSE(std::filesystem::exists(untouched));
    EXPECT_EQ(static_cast<std::string>(same["LATTICE"]), "square lattice");
    EXPECT_THROW(lat::write_lattice_library(resolution, std::filesystem::path("/nonexistent/dir/x.xml")), std::runtime_error);
}

TEST(LatticeXmlBridge, CatalogModelOnAKitLatticeReproducesTheKitGroundState) {
    // The typed honeycomb cell needs the custom route; the catalog model
    // addresses its three bond types.
    auto const honeycomb = lat::honeycomb(2, 2);
    alps::Parameters p;
    p["MODEL"] = "spin";
    p["local_S"] = 0.5;
    p["J0"] = 1.0;
    p["J1"] = 1.0;
    p["J2"] = 1.0;
    mdl::CatalogModel<> model(honeycomb, p);
    EXPECT_EQ(static_cast<int>(model.graph().num_sites()), honeycomb.num_sites());
    EXPECT_EQ(static_cast<int>(model.graph().num_bonds()), honeycomb.num_bonds());
    EXPECT_EQ(static_cast<int>(model.graph().site_type(1)), 1);
    EXPECT_EQ(static_cast<std::string>(model.parameters()["LATTICE"]), lat::custom_lattice_name(honeycomb));
    EXPECT_NEAR(catalog_e0(model), kit_e0(mdl::heisenberg(honeycomb, 1.0)), 1e-9);
    // Distinct leg and rung couplings on an open ladder.
    auto const ladder = lat::ladder(3, 1.0, Boundary::open);
    alps::Parameters q;
    q["MODEL"] = "spin";
    q["local_S"] = 0.5;
    q["J0"] = 1.0;
    q["J1"] = 0.4;
    mdl::CatalogModel<> rungs(ladder, q);
    auto const expected = mdl::heisenberg(ladder, 1.0, 0.0, 0.5, 0);
    mdl::Hamiltonian reference = expected;
    auto const site = reference.site_for(0);
    reference.add_bond(0.4, [site](int i, int j) {
        using mdl::TermSum;
        return TermSum::located(site, "Sz", i) * TermSum::located(site, "Sz", j)
             + std::complex<double>(0.5) * (TermSum::located(site, "Sp", i) * TermSum::located(site, "Sm", j)
                                            + TermSum::located(site, "Sm", i) * TermSum::located(site, "Sp", j));
    }, 1);
    EXPECT_NEAR(catalog_e0(rungs), kit_e0(reference), 1e-9);
}
