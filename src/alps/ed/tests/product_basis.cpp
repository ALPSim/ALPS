// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/tests/models/product_basis.cpp.
// Adaptation: <alps/ed/...> include names.

#include <alps/ed/product_basis.hpp>
#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <stdexcept>

using alps::ed::BasisConstraint;
using alps::ed::enumerate_product_basis;

TEST(ProductBasis, FullOrderAndVacuum) {
    std::vector<int> const dims{2, 3};
    auto const full = enumerate_product_basis(dims);
    EXPECT_EQ(full.dimension, 6);
    EXPECT_EQ(full.states, (std::vector<std::int32_t>{0,0, 0,1, 0,2, 1,0, 1,1, 1,2}));
    auto const vacuum = enumerate_product_basis({});
    EXPECT_EQ(vacuum.dimension, 1);
    EXPECT_TRUE(vacuum.states.empty());
    std::vector<BasisConstraint> constraints{{{}, 1.}};
    EXPECT_EQ(enumerate_product_basis({}, constraints).dimension, 0);
    constraints[0].target = 0.;
    EXPECT_EQ(enumerate_product_basis({}, constraints).dimension, 1);
}

TEST(ProductBasis, MultipleChargesAndStrictTolerance) {
    std::vector<int> const dims{3, 2, 2};
    std::vector<BasisConstraint> const constraints{
        {{{0., 1., 2.}, {0., 1.}, {0., 1.}}, 2.},
        {{{-.5, .5, 1.5}, {-.5, .5}, {-.5, .5}}, .5}};
    auto const basis = enumerate_product_basis(dims, constraints);
    EXPECT_EQ(basis.states, (std::vector<std::int32_t>{0,1,1, 1,0,1, 1,1,0, 2,0,0}));
    EXPECT_EQ(basis.dimension, 4);
    std::vector<int> const one{4};
    std::vector<BasisConstraint> const boundary{
        {{{0., std::nextafter(1e-9, 0.), 1e-9, std::nextafter(1e-9, 1.)}}, 0.}};
    EXPECT_EQ(enumerate_product_basis(one, boundary).states,
              (std::vector<std::int32_t>{0, 1}));
    std::vector<int> const scalar_sites{1, 1, 1};
    std::vector<BasisConstraint> const rounding{{{{1e16}, {-1e16}, {1.}}, 1.}};
    EXPECT_EQ(enumerate_product_basis(scalar_sites, rounding).dimension, 1);
}

TEST(ProductBasis, DeepProjectedSpaceDoesNotEnumerateFullProductOrUseRecursion) {
    std::vector<int> const dims(2000, 2);
    std::vector<BasisConstraint> const constraints{
        {std::vector<std::vector<double>>(dims.size(), {0., 1.}), 0.}};
    auto const basis = enumerate_product_basis(dims, constraints);
    EXPECT_EQ(basis.dimension, 1);
    EXPECT_EQ(basis.states, std::vector<std::int32_t>(dims.size(), 0));
}

TEST(ProductBasis, LargeLocalIndicesAndCheckedFullDimension) {
    std::vector<int> const dims{40000};
    auto const full = enumerate_product_basis(dims);
    EXPECT_EQ(full.dimension, 40000);
    EXPECT_EQ(full.states.front(), 0);
    EXPECT_EQ(full.states.back(), 39999);
    std::vector<int> const enormous(64, 2);
    EXPECT_THROW(enumerate_product_basis(enormous), std::invalid_argument);
    EXPECT_THROW(enumerate_product_basis(dims, {}, 39999), std::invalid_argument);
    EXPECT_THROW(enumerate_product_basis({}, {}, 0), std::invalid_argument);
}

TEST(ProductBasis, RejectsMalformedChargeArrays) {
    std::vector<int> const dims{2};
    std::vector<int> const invalid{0};
    EXPECT_THROW(enumerate_product_basis(invalid), std::invalid_argument);
    std::vector<BasisConstraint> constraints{{{}, 0.}};
    EXPECT_THROW(enumerate_product_basis(dims, constraints), std::invalid_argument);
    constraints[0].values = {{0.}};
    EXPECT_THROW(enumerate_product_basis(dims, constraints), std::invalid_argument);
    constraints[0].values = {{0., std::numeric_limits<double>::quiet_NaN()}};
    EXPECT_THROW(enumerate_product_basis(dims, constraints), std::invalid_argument);
    constraints[0].values = {{0., 1.}};
    constraints[0].target = std::numeric_limits<double>::infinity();
    EXPECT_THROW(enumerate_product_basis(dims, constraints), std::invalid_argument);
}
