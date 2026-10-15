// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// FiniteGraph and LatticeGeometry are owned inputs without a unit-cell
// tiler: these cases cover their validation contract. Construction of a
// catalog graph from them is a consumer of the legacy lattice library and is
// tested there.

#include <alps/models/finite_graph.hpp>
#include <alps/models/geometry.hpp>

#include <gtest/gtest.h>

#include <stdexcept>

namespace {

alps::FiniteGraph triangle() {
    alps::FiniteGraph g;
    g.name = "triangle";
    g.dimension = 2;
    g.vertices = {{0, {0., 0.}}, {0, {1., 0.}}, {1, {.5, .5}}};
    g.edges = {{0, 1, 0, {1., 0.}}, {1, 2, 1, {-.5, .5}}, {2, 0, 1, {-.5, -.5}}};
    return g;
}

} // namespace

TEST(FiniteGraph, ValidatesAWellFormedGraph) {
    EXPECT_NO_THROW(triangle().validate());
    auto plain = triangle();
    plain.dimension = 0;
    for (auto& v : plain.vertices) v.coordinate.clear();
    for (auto& e : plain.edges) e.displacement.clear();
    EXPECT_NO_THROW(plain.validate());
}

TEST(FiniteGraph, RejectsMalformedInputs) {
    auto empty = triangle();
    empty.vertices.clear();
    empty.edges.clear();
    EXPECT_THROW(empty.validate(), std::invalid_argument);

    auto dangling = triangle();
    dangling.edges[0].target = 7;
    EXPECT_THROW(dangling.validate(), std::invalid_argument);

    auto mismatched = triangle();
    mismatched.vertices[1].coordinate = {1.};
    EXPECT_THROW(mismatched.validate(), std::invalid_argument);
}

TEST(FiniteGraph, ComparesByValue) {
    EXPECT_EQ(triangle(), triangle());
    auto other = triangle();
    other.edges[1].type = 3;
    EXPECT_NE(triangle(), other);
}

TEST(LatticeGeometry, ValidatesExtentBoundaryAndCell) {
    alps::LatticeGeometry geometry;
    geometry.extent = {2, 2};
    geometry.boundary = {"periodic", "open"};
    geometry.vertices = {{0, {0., 0.}}};
    geometry.edges = {{0, 0, 0, {0, 0}, {1, 0}}, {0, 0, 0, {0, 0}, {0, 1}}};
    EXPECT_NO_THROW(geometry.validate());

    auto no_cell = geometry;
    no_cell.vertices.clear();
    EXPECT_THROW(no_cell.validate(), std::invalid_argument);

    auto bad_offset = geometry;
    bad_offset.edges[0].target_offset = {1};
    EXPECT_THROW(bad_offset.validate(), std::invalid_argument);
}
