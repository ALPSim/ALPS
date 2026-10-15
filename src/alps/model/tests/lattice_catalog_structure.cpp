// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Structural comparison of the kit lattice factories (ALPS::models) with the
// catalog graphs of the same name (lattices.xml through graph_helper): site
// and bond counts, the bond multiset with its types, site types and vertex
// coordinates. Known divergences are asserted explicitly so that a change on
// either side is noticed.
#include <alps/lattice.h>
#include <alps/models/lattice.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace lat = alps::lattice;

namespace {
using Edges = std::multiset<std::tuple<int, int, int>>; // (min site, max site, bond type)

Edges kit_edges(lat::Lattice const& lattice) {
    Edges edges;
    for (auto const& bond : lattice.bonds())
        edges.emplace(std::min(bond.source, bond.target), std::max(bond.source, bond.target), bond.bond_type);
    return edges;
}
Edges catalog_edges(alps::graph_helper<> const& graph) {
    Edges edges;
    for (auto it = graph.bonds().first; it != graph.bonds().second; ++it) {
        int const s = static_cast<int>(graph.source(*it)), t = static_cast<int>(graph.target(*it));
        edges.emplace(std::min(s, t), std::max(s, t), static_cast<int>(graph.bond_type(*it)));
    }
    return edges;
}
std::vector<int> catalog_site_types(alps::graph_helper<> const& graph) {
    std::vector<int> types;
    for (auto it = graph.sites().first; it != graph.sites().second; ++it)
        types.push_back(static_cast<int>(graph.site_type(*it)));
    return types;
}
std::vector<std::vector<double>> catalog_coordinates(alps::graph_helper<> const& graph) {
    std::vector<std::vector<double>> coordinates;
    for (auto it = graph.sites().first; it != graph.sites().second; ++it) {
        auto const& c = boost::get(alps::coordinate_t(), graph.graph(), *it);
        coordinates.emplace_back(c.begin(), c.end());
    }
    return coordinates;
}
alps::Parameters parameters(std::string const& name, lat::Lattice const& lattice) {
    alps::Parameters p;
    p["LATTICE"] = name;
    p["L"] = lattice.extent().at(0);
    if (lattice.extent().size() > 1) p["W"] = lattice.extent().at(1);
    if (lattice.extent().size() > 2) p["H"] = lattice.extent().at(2);
    return p;
}

struct Expectation {
    bool site_types_match = true;  // false: the catalog leaves every site untyped
    bool coordinates_match = true; // false: documented coordinate divergence
};

// `name` is the catalog LATTICEGRAPH; the kit's alps_name() names the
// periodic graph for every boundary except the chains, so the open variants
// are spelled out here.
void compare(lat::Lattice const& lattice, std::string const& name, Expectation const& expect = {}) {
    SCOPED_TRACE(name);
    alps::graph_helper<> const graph(parameters(name, lattice));
    EXPECT_EQ(static_cast<int>(graph.num_sites()), lattice.num_sites());
    EXPECT_EQ(static_cast<int>(graph.num_bonds()), lattice.num_bonds());
    EXPECT_EQ(catalog_edges(graph), kit_edges(lattice));
    auto const types = catalog_site_types(graph), kit_types = lattice.site_types();
    if (expect.site_types_match) {
        EXPECT_EQ(types, kit_types);
    } else {
        EXPECT_TRUE(std::all_of(types.begin(), types.end(), [](int t) { return t == 0; }));
        EXPECT_GT(*std::max_element(kit_types.begin(), kit_types.end()), 0);
    }
    auto const coordinates = catalog_coordinates(graph);
    ASSERT_EQ(coordinates.size(), lattice.sites().size());
    double distance = 0.;
    for (std::size_t site = 0; site < coordinates.size(); ++site) {
        auto const& kit = lattice.sites()[site].coordinate;
        ASSERT_EQ(coordinates[site].size(), kit.size());
        for (std::size_t d = 0; d < kit.size(); ++d)
            distance = std::max(distance, std::abs(coordinates[site][d] - kit[d]));
    }
    if (expect.coordinates_match) EXPECT_LE(distance, 1e-12);
    else EXPECT_GT(distance, 1e-12);
}
} // namespace

TEST(LatticeCatalogStructure, SingleVertexCellsMatchTheCatalogExactly) {
    using B = lat::Boundary;
    compare(lat::chain(5), "chain lattice");
    compare(lat::chain(5, 1.0, B::open), "open chain lattice");
    compare(lat::nnn_chain(5), "nnn chain lattice");
    compare(lat::nnn_chain(5, 1.0, B::open), "nnn open chain lattice");
    compare(lat::ladder(4), "ladder");
    compare(lat::ladder(4, 1.0, B::open), "open ladder");
    compare(lat::square(3, 4), "square lattice");
    compare(lat::square(3, 4, 1.0, B::open), "open square lattice");
    compare(lat::frustrated_square(3, 3), "frustrated square lattice");
    compare(lat::triangular(3, 4), "triangular lattice");
    compare(lat::simple_cubic(2, 3, 2), "simple cubic lattice");
    compare(lat::simple_cubic(2, 3, 2, 1.0, B::open), "open simple cubic lattice");
}

TEST(LatticeCatalogStructure, KitNamesTheOpenVariantsOnlyForChains) {
    using B = lat::Boundary;
    EXPECT_EQ(lat::chain(4, 1.0, B::open).alps_name(), "open chain lattice");
    EXPECT_EQ(lat::nnn_chain(4, 1.0, B::open).alps_name(), "nnn open chain lattice");
    EXPECT_EQ(lat::ladder(4, 1.0, B::open).alps_name(), "ladder");
    EXPECT_EQ(lat::square(3, 3, 1.0, B::open).alps_name(), "square lattice");
    EXPECT_EQ(lat::simple_cubic(2, 2, 2, 1.0, B::open).alps_name(), "simple cubic lattice");
}

TEST(LatticeCatalogStructure, MultiVertexCellsDivergeInSiteTypes) {
    // The kit tags the sublattices (0/1 honeycomb, 0/1/2 kagome); the catalog
    // leaves every vertex of these unit cells untyped.
    Expectation untyped;
    untyped.site_types_match = false;
    compare(lat::honeycomb(2, 3), "honeycomb lattice", untyped);
}

TEST(LatticeCatalogStructure, KagomeNeedsTheCatalogCasing) {
    auto const kagome = lat::kagome(2, 2);
    EXPECT_EQ(kagome.alps_name(), "kagome lattice");
    EXPECT_THROW(alps::graph_helper<>(parameters(kagome.alps_name(), kagome)), std::runtime_error);
    alps::graph_helper<> const graph(parameters("Kagome lattice", kagome));
    EXPECT_EQ(static_cast<int>(graph.num_sites()), kagome.num_sites());
    EXPECT_EQ(catalog_edges(graph), kit_edges(kagome));
    auto const types = catalog_site_types(graph);
    EXPECT_TRUE(std::all_of(types.begin(), types.end(), [](int t) { return t == 0; }));
}

TEST(LatticeCatalogStructure, FactoriesWithoutACatalogCounterpart) {
    for (auto const& lattice : {lat::rectangular(3, 2, 1.0, 2.0), lat::body_centered_cubic(2),
                                lat::face_centered_cubic(2), lat::tetragonal(2), lat::orthorhombic(2),
                                lat::hexagonal(2)}) {
        SCOPED_TRACE(lattice.alps_name());
        EXPECT_THROW(alps::graph_helper<>(parameters(lattice.alps_name(), lattice)), std::runtime_error);
    }
}
