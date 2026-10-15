// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/tests/models/streamed_local_operator.cpp.
// Adaptation: C++17 (no std::span), <alps/ed/...> include names.
#include <alps/ed/sparse_from_local.hpp>
#include <gtest/gtest.h>

#include <complex>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace {
using Complex = std::complex<double>;
using Operator = alps::models::CompiledOperator;
using Coo = alps::ed::CooFromLocal;
using Visitor = alps::ed::MatrixElementVisitor;

Visitor collect(Coo& out) {
    return [&out](std::int64_t row, std::int64_t col, Complex value) {
        out.rows.push_back(row);
        out.cols.push_back(col);
        out.data.push_back(value);
    };
}

void same(Coo const& actual, Coo const& expected) {
    EXPECT_EQ(actual.rows, expected.rows);
    EXPECT_EQ(actual.cols, expected.cols);
    EXPECT_EQ(actual.data, expected.data);
    EXPECT_EQ(actual.is_complex, expected.is_complex);
}
} // namespace

TEST(StreamedLocalOperator, RetainsKnownContributionOrderProjectionAndCancellation) {
    Operator op{{2, 2}, {
        {{0., 2.}, {}, {}},
        {{1., .5}, {{0, {2., 0., 0., -1.}}, {1, {0., 3., 4., 0.}}}, {0}},
        {{-1., -.5}, {{0, {2., 0., 0., -1.}}, {1, {0., 3., 4., 0.}}}, {0}},
        {1., {{0, {0., 0., 5., 0.}}}, {}},
        {.25, {{0, {1., 2., 3., 4.}}}, {1}}}};
    std::vector<std::int64_t> const states{1, 0, 0, 1, 1, 1};
    std::vector<std::uint64_t> const codes{2, 1, 3};
    std::vector<std::vector<double>> const parity{{1., -1.}, {0., -1.}};
    Coo const expected{{0, 1, 2, 2, 0, 2, 0, 2, 1, 2, 1, 2},
                       {0, 1, 2, 0, 2, 0, 2, 1, 1, 1, 2, 2},
                       {{0., 2.}, {0., 2.}, {0., 2.}, {4., 2.}, {3., 1.5},
                        {-4., -2.}, {-3., -1.5}, 5., -.25, -.75, -.5, -1.}, true};
    Coo explicit_rows, encoded;
    explicit_rows.is_complex = alps::ed::visit_local_operator(op, states, 3, collect(explicit_rows), parity);
    encoded.is_complex = alps::ed::visit_local_operator(op, codes, collect(encoded), parity);
    same(explicit_rows, expected);
    same(encoded, expected);
    same(alps::ed::build_coo_from_local(op, states, 3, parity), expected);
}

TEST(StreamedLocalOperator, SortedMixedRadixCodesUseSiteZeroMostSignificant) {
    Operator op{{2, 3}, {{2., {{1, {0., 1., 0., 0., 0., 3., 4., 0., 0.}}}, {0}}}};
    std::vector<std::uint64_t> const codes{0, 1, 2, 3, 4, 5};
    std::vector<std::int64_t> const states{0, 0, 0, 1, 0, 2, 1, 0, 1, 1, 1, 2};
    std::vector<std::vector<double>> const parity{{1., -1.}, {1., 1., 1.}};
    Coo encoded;
    encoded.is_complex = alps::ed::visit_local_operator(op, codes, collect(encoded), parity);
    Coo const expected{{2, 0, 1, 5, 3, 4}, {0, 1, 2, 3, 4, 5}, {8., 2., 6., -8., -2., -6.}, false};
    same(encoded, expected);
    same(alps::ed::build_coo_from_local(op, states, 6, parity), expected);
}

TEST(StreamedLocalOperator, UnsignedCodesRetainTheHighBitWithoutSignedOverflow) {
    Operator op;
    op.site_dims.assign(63, 2);
    op.site_dims[0] = 3;
    op.terms = {{1., {{0, {0., 0., 1., 0., 0., 0., 0., 0., 0.}}}, {}}};
    std::vector<std::uint64_t> const codes{std::uint64_t{1} << 63, 0};
    Coo encoded;
    encoded.is_complex = alps::ed::visit_local_operator(op, codes, collect(encoded));
    same(encoded, Coo{{1}, {0}, {1.}, false});
}

TEST(StreamedLocalOperator, ExplicitRowsSupportProductsThatCannotBeEncoded) {
    Operator op;
    op.site_dims.assign(70, 2);
    op.terms = {{{-.25, -.5}, {{0, {0., 1., 0., 0.}}, {69, {0., 0., 1., 0.}}}, {}}};
    for (int site = 0; site < 69; ++site) op.terms[0].jordan_wigner.push_back(site);
    std::vector<std::int64_t> states(140, 0);
    states[0] = states[35] = states[105] = states[139] = 1;
    std::vector<std::vector<double>> const parity(70, {1., -1.});
    Coo explicit_rows;
    explicit_rows.is_complex = alps::ed::visit_local_operator(op, states, 2, collect(explicit_rows), parity);
    same(explicit_rows, Coo{{1}, {0}, {{-.25, -.5}}, true});
    std::vector<std::uint64_t> const codes{0};
    EXPECT_THROW(alps::ed::visit_local_operator(op, codes, collect(explicit_rows), parity), std::invalid_argument);
    op.site_dims.assign(64, 2); // Product 2^64 is too large, even though its maximum code fits.
    op.terms.clear();
    EXPECT_THROW(alps::ed::visit_local_operator(op, codes, collect(explicit_rows)), std::invalid_argument);
}

TEST(StreamedLocalOperator, ValidatesBothRepresentationsBeforeAnyCallback) {
    Operator const valid{{2}, {{1., {{0, {0., 1., 1., 0.}}}, {}}}};
    std::size_t callbacks = 0;
    Visitor const visitor = [&](auto, auto, auto) { ++callbacks; };
    for (auto const& codes : {std::vector<std::uint64_t>{0, 0}, {1, 0, 1}, {2}})
        EXPECT_THROW(alps::ed::visit_local_operator(valid, codes, visitor), std::invalid_argument);
    std::vector<std::uint64_t> const codes{0, 1};
    std::vector<std::int64_t> const states{0, 1};
    auto bad = valid;
    bad.site_dims[0] = 0;
    EXPECT_THROW(alps::ed::visit_local_operator(bad, codes, visitor), std::invalid_argument);
    EXPECT_THROW(alps::ed::visit_local_operator(bad, states, 2, visitor), std::invalid_argument);
    for (auto const& parity : {std::vector<std::vector<double>>{{1.}}, {{1., -1.}, {1., -1.}}}) {
        EXPECT_THROW(alps::ed::visit_local_operator(valid, codes, visitor, parity), std::invalid_argument);
        EXPECT_THROW(alps::ed::visit_local_operator(valid, states, 2, visitor, parity), std::invalid_argument);
    }
    bad = valid;
    bad.terms.push_back({1., {{0, {1.}}}, {}}); // A valid first term must not escape first.
    EXPECT_THROW(alps::ed::visit_local_operator(bad, codes, visitor), std::invalid_argument);
    EXPECT_THROW(alps::ed::visit_local_operator(bad, states, 2, visitor), std::invalid_argument);
    bad = valid;
    bad.terms[0].jordan_wigner = {1};
    EXPECT_THROW(alps::ed::visit_local_operator(bad, codes, visitor), std::invalid_argument);
    EXPECT_THROW(alps::ed::visit_local_operator(valid, states, 3, visitor), std::invalid_argument);
    EXPECT_EQ(callbacks, 0);
}

TEST(StreamedLocalOperator, PropagatesCallbackExceptionsAndRejectsEmptyVisitors) {
    Operator const op{{2}, {{1., {}, {}}, {2., {}, {}}}};
    std::vector<std::uint64_t> const codes{0, 1};
    std::vector<std::int64_t> const states{0, 1};
    struct CallbackError {};
    for (bool encoded : {false, true}) {
        std::size_t callbacks = 0;
        Visitor const visitor = [&](auto, auto, auto) { ++callbacks; throw CallbackError{}; };
        if (encoded) EXPECT_THROW(alps::ed::visit_local_operator(op, codes, visitor), CallbackError);
        else EXPECT_THROW(alps::ed::visit_local_operator(op, states, 2, visitor), CallbackError);
        EXPECT_EQ(callbacks, 1);
    }
    EXPECT_THROW(alps::ed::visit_local_operator(op, codes, {}), std::invalid_argument);
    EXPECT_THROW(alps::ed::visit_local_operator(op, states, 2, {}), std::invalid_argument);
    EXPECT_THROW(alps::ed::visit_local_operator(op, std::vector<std::uint64_t>{}, {}), std::invalid_argument);
}

TEST(StreamedLocalOperator, ComplexFlagMatchesEmittedContributionsIncludingCoefficientPhase) {
    std::vector<std::uint64_t> const codes{0};
    Operator op{{1}, {{{0., 1.}, {{0, {{0., -1.}}}}, {}}}};
    Coo out;
    out.is_complex = alps::ed::visit_local_operator(op, codes, collect(out));
    same(out, Coo{{0}, {0}, {1.}, true}); // Real value from complex coefficient still marks the stream.
    op.terms[0] = {1., {{0, {{0., 1e-100}}}}, {}};
    out = {};
    out.is_complex = alps::ed::visit_local_operator(op, codes, collect(out));
    same(out, Coo{{0}, {0}, {{0., 1e-100}}, true});
    op.site_dims = {2};
    op.terms[0] = {{0., 1.}, {{0, {0., 0., 1., 0.}}}, {}};
    out = {};
    out.is_complex = alps::ed::visit_local_operator(op, codes, collect(out));
    same(out, Coo{}); // The complex transition leaves this selected sector.
}

TEST(StreamedLocalOperator, ZeroSiteVacuumAndEmptySectorRemainDistinct) {
    Operator const op{{}, {{{0., 2.}, {}, {}}}};
    std::vector<std::uint64_t> const vacuum{0};
    Coo out;
    out.is_complex = alps::ed::visit_local_operator(op, vacuum, collect(out));
    same(out, Coo{{0}, {0}, {{0., 2.}}, true});
    out = {};
    out.is_complex = alps::ed::visit_local_operator(op, std::vector<std::uint64_t>{}, collect(out));
    same(out, Coo{});
    EXPECT_THROW(alps::ed::visit_local_operator(op, std::vector<std::uint64_t>{1}, collect(out)), std::invalid_argument);
    EXPECT_THROW(alps::ed::visit_local_operator(op, std::vector<std::uint64_t>{0, 0}, collect(out)), std::invalid_argument);
    EXPECT_THROW(alps::ed::visit_local_operator(op, std::vector<std::int64_t>{}, 2, collect(out)), std::invalid_argument);
}
