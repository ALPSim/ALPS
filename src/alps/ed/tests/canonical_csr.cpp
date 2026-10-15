// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/tests/canonical_csr.cpp.
// Adaptation: <alps/ed/...> include names.
#include <alps/ed/sparse_matrix.hpp>
#include <alps/models/builtin.hpp>
#include <gtest/gtest.h>

#include <cstdint>
#include <limits>
#include <stdexcept>

namespace {
using alps::ed::CooFromLocal;
using alps::ed::canonical_csr;
using C = std::complex<double>;
}
TEST(CanonicalCsr, RectangularReductionAndOwnership) {
    CooFromLocal coo{{1,0,1,0,1}, {2,1,2,0,0}, {{2,1},3.,{-2,-1},4.,5.}, true};
    auto result = canonical_csr<C>(coo, 2, 3);
    coo.data.assign(5, 99.);
    EXPECT_EQ(result.rows, 2); EXPECT_EQ(result.columns, 3);
    EXPECT_EQ(result.offsets, (std::vector<std::int64_t>{0,2,3}));
    EXPECT_EQ(result.indices, (std::vector<std::int64_t>{0,1,0}));
    EXPECT_EQ(result.values, (std::vector<C>{4.,3.,5.}));
    auto copy = result;
    result.values[0] = 99.;
    EXPECT_EQ(copy.values[0], C(4.));
    auto empty = canonical_csr<double>({}, 0, 7);
    EXPECT_EQ(empty.offsets, (std::vector<std::int64_t>{0}));
    EXPECT_TRUE(empty.values.empty());
}
TEST(CanonicalCsr, ConversionPrecedesReductionAndRealLossIsRejected) {
    CooFromLocal coo{{0,0,0},{0,0,0},{16777216.,1.,-16777216.}, false};
    EXPECT_EQ(canonical_csr<double>(coo,1,1).values, (std::vector<double>{1.}));
    EXPECT_TRUE(canonical_csr<float>(coo,1,1).values.empty());
    EXPECT_TRUE(canonical_csr<std::complex<float>>(coo,1,1).values.empty());
    coo = {{0,0},{0,0},{{0.,1.},{0.,-1.}},true};
    EXPECT_THROW(canonical_csr<double>(coo,1,1), std::invalid_argument);
    EXPECT_TRUE(canonical_csr<C>(coo,1,1).values.empty());
    coo = {{0},{0},{{1.,1e-100}},true};
    EXPECT_THROW(canonical_csr<float>(coo,1,1), std::invalid_argument);
}
TEST(CanonicalCsr, DefinedIntegerAndBooleanReduction) {
    CooFromLocal coo{{0,0},{0,0},{127.,1.}, false};
    EXPECT_EQ(canonical_csr<std::int8_t>(coo,1,1).values, (std::vector<std::int8_t>{-128}));
    EXPECT_EQ(canonical_csr<bool>(coo,1,1).values, (std::vector<bool>{true}));
    coo.data = {-1.,1.};
    EXPECT_TRUE(canonical_csr<std::uint64_t>(coo,1,1).values.empty());
    coo.data = {-1.9,0.1};
    EXPECT_EQ(canonical_csr<std::int64_t>(coo,1,1).values, (std::vector<std::int64_t>{-1}));
}
TEST(CanonicalCsr, RejectsMalformedAndNonfiniteContributions) {
    EXPECT_THROW(canonical_csr<double>({{0},{},{1.}, false},1,1), std::invalid_argument);
    EXPECT_THROW(canonical_csr<double>({{-1},{0},{1.}, false},1,1), std::invalid_argument);
    EXPECT_THROW(canonical_csr<double>({{0},{1},{1.}, false},1,1), std::invalid_argument);
    EXPECT_THROW(canonical_csr<double>({{0},{0},{std::numeric_limits<double>::infinity()}, false},1,1), std::invalid_argument);
    EXPECT_THROW(canonical_csr<float>({{0},{0},{1e100}, false},1,1), std::invalid_argument);
    EXPECT_THROW(canonical_csr<double>({{0,0},{0,0},{1e308,1e308}, false},1,1), std::invalid_argument);
    EXPECT_THROW(canonical_csr<double>({},std::numeric_limits<std::size_t>::max(),1), std::length_error);
}
TEST(CanonicalCsr, ArbitraryOperatorPreservesCallerRows) {
    namespace mdl = alps::models;
    auto spin = mdl::spin();
    auto op = mdl::compile_operator(mdl::TermSum::located(spin,"Sp",0), {spin,spin});
    std::vector<std::int64_t> states{1,0,0,0,0,1};
    auto result = alps::ed::build_csr_from_local<double>(op, states, 3);
    EXPECT_EQ(result.offsets, (std::vector<std::int64_t>{0,0,1,1}));
    EXPECT_EQ(result.indices, (std::vector<std::int64_t>{0}));
    EXPECT_EQ(result.values, (std::vector<double>{1.}));
}
