// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/tests/lattice/translation_generators.cpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports;
// TOML emitter check replaced by a basis-vector check.
#include <alps/models/lattice.hpp>
#include <gtest/gtest.h>
#include <stdexcept>

TEST(TranslationGenerators, SourceToDestinationAndOpenDirections) {
    namespace lat = alps::lattice;
    EXPECT_EQ(lat::chain(4).translation_generators(),
              (std::vector<std::vector<std::size_t>>{{1, 2, 3, 0}}));
    EXPECT_TRUE(lat::chain(4, 1., lat::Boundary::open).translation_generators().empty());
    EXPECT_EQ(lat::chain(1).translation_generators(),
              (std::vector<std::vector<std::size_t>>{{0}}));
    EXPECT_EQ(lat::square(2, 3).translation_generators(),
              (std::vector<std::vector<std::size_t>>{{3,4,5,0,1,2}, {1,2,0,4,5,3}}));
}

TEST(TranslationGenerators, CoincidentVerticesKeepTheirIdentity) {
    namespace lat = alps::lattice;
    lat::UnitCell cell(2);
    cell.add_vertex({0., 0.}, 0);
    cell.add_vertex({0., 0.}, 1);
    auto cylinder = cell.tile({3, 2}, {{1., .2}, {0., 2.}},
                              {lat::Boundary::periodic, lat::Boundary::open});
    EXPECT_EQ(cylinder.periodic_directions(), (std::vector<int>{0}));
    EXPECT_EQ(cylinder.translation_generators(),
              (std::vector<std::vector<std::size_t>>{{4,5,6,7,8,9,10,11,0,1,2,3}}));
}

TEST(TranslationGenerators, EmbeddedCoordinatesAndWrappedBondVectors) {
    namespace lat = alps::lattice;
    lat::UnitCell cell(2);
    cell.add_vertex({.25, .5}, 7);
    cell.add_edge(0, 0, 3, {}, {1, 0});
    auto lattice = cell.tile({2, 2}, {{2., 0., 1.}, {0., 4., 2.}},
                             {lat::Boundary::periodic, lat::Boundary::open});
    std::vector<std::vector<double>> const coordinates{
        {.5, 2., 1.25}, {.5, 6., 3.25}, {2.5, 2., 2.25}, {2.5, 6., 4.25}};
    ASSERT_EQ(lattice.num_sites(), 4);
    ASSERT_EQ(lattice.num_bonds(), 4);
    for (std::size_t i = 0; i < coordinates.size(); ++i) {
        EXPECT_EQ(lattice.sites()[i].coordinate, coordinates[i]);
        EXPECT_EQ(lattice.sites()[i].site_type, 7);
        EXPECT_EQ(lattice.bonds()[i].vector,
                  (i < 2 ? std::vector<double>{2., 0., 1.} : std::vector<double>{-2., 0., -1.}));
        EXPECT_EQ(lattice.bonds()[i].wraps, i >= 2);
    }
    EXPECT_EQ(lattice.dimension(), 2);
    EXPECT_EQ(lattice.translation_generators(),
              (std::vector<std::vector<std::size_t>>{{2, 3, 0, 1}}));
    EXPECT_EQ(lattice.basis_vectors(),
              (std::vector<std::vector<double>>{{2., 0., 1.}, {0., 4., 2.}}));
}

TEST(TranslationGenerators, RejectsInconsistentConstructionShapes) {
    namespace lat = alps::lattice;
    lat::UnitCell cell(2);
    cell.add_vertex({0., 0.});
    auto const boundary = std::vector<lat::Boundary>(2, lat::Boundary::periodic);
    EXPECT_THROW(lat::Lattice::build(cell, {{1., 0.}, {0., 1.}}, {2}, boundary),
                 std::invalid_argument);
    EXPECT_THROW(lat::Lattice::build(cell, {{1., 0.}, {0., 1.}}, {2, 2}, {}),
                 std::invalid_argument);
    EXPECT_THROW(lat::Lattice::build(cell, {{1., 0.}}, {2, 2}, boundary),
                 std::invalid_argument);
    EXPECT_THROW(lat::Lattice::build(cell, {{1., 0.}, {0.}}, {2, 2}, boundary),
                 std::invalid_argument);
    EXPECT_THROW(lat::Lattice::build(cell, {{}, {}}, {2, 2}, boundary),
                 std::invalid_argument);
}
