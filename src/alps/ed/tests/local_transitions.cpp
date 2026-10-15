// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/tests/models/local_transitions.cpp.
// Adaptation: C++17 (no std::span or defaulted comparisons), <alps/ed/...> include names.
#include <alps/ed/sparse_from_local.hpp>
#include <gtest/gtest.h>

#include <complex>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
using Complex = std::complex<double>;
using Operator = alps::models::CompiledOperator;
struct Transition {
    std::int64_t source;
    std::vector<std::int64_t> target;
    Complex amplitude;
    bool operator==(Transition const& o) const {
        return source == o.source && target == o.target && amplitude == o.amplitude;
    }
};
auto collect(std::vector<Transition>& transitions) {
    return [&transitions](std::int64_t source, std::vector<std::int64_t> const& target, Complex amplitude) {
        transitions.push_back({source, target, amplitude});
    };
}
} // namespace

TEST(LocalTransitions, EmitsOrderedTargetsOutsideSourceBatchWithParityAndBranching) {
    Operator const op{{2, 3}, {
        {2., {{1, {0., 1., 0., 0., 0., 3., 4., 0., 0.}}}, {0}},
        {.5, {{0, {1., 2., 3., 4.}}, {1, {1., 0., 0., 0., 2., 0., 0., 0., 3.}}}, {1}}}};
    std::vector<std::int64_t> const sources{1, 0, 0, 2};
    std::vector<std::vector<double>> const parity{{1., -1.}, {1., -1., 1.}};
    std::vector<Transition> actual;
    EXPECT_FALSE(alps::ed::visit_local_transitions(op, sources, 2, collect(actual), parity));
    std::vector<Transition> const expected{
        {0, {1, 2}, -8.}, {1, {0, 1}, 6.},
        {0, {0, 0}, 1.}, {0, {1, 0}, 2.},
        {1, {0, 2}, 1.5}, {1, {1, 2}, 4.5}};
    EXPECT_EQ(actual, expected);
}

TEST(LocalTransitions, RetainsDuplicateSourcesAndComplexTermsWithoutAggregationOrCutoff) {
    Operator const op{{1}, {
        {{0., 2.}, {{0, {{0., -1.}}}}, {}},
        {-2., {}, {}}, {{0., 1e-100}, {}, {}}, {0., {}, {}}}};
    std::vector<std::int64_t> const sources{0, 0};
    std::vector<Transition> actual;
    EXPECT_TRUE(alps::ed::visit_local_transitions(op, sources, 2, collect(actual)));
    std::vector<Transition> const expected{
        {0, {0}, 2.}, {1, {0}, 2.}, {0, {0}, -2.}, {1, {0}, -2.},
        {0, {0}, {0., 1e-100}}, {1, {0}, {0., 1e-100}}};
    EXPECT_EQ(actual, expected);
    Operator const real_amplitude{{1}, {{{0., 2.}, {{0, {{0., -1.}}}}, {}}}};
    actual.clear();
    EXPECT_TRUE(alps::ed::visit_local_transitions(real_amplitude, sources, 2, collect(actual)));
    EXPECT_EQ(actual.front().amplitude, 2.);
}

TEST(LocalTransitions, SupportsUnencodableProductWithoutEnumeratingIt) {
    Operator op;
    op.site_dims.assign(70, 2);
    op.terms = {{{-.25, -.5}, {{0, {0., 1., 0., 0.}}, {69, {0., 0., 1., 0.}}}, {35}}};
    std::vector<std::int64_t> source(70);
    source[0] = source[35] = 1;
    std::vector<std::vector<double>> const parity(70, {1., -1.});
    std::vector<Transition> actual;
    EXPECT_TRUE(alps::ed::visit_local_transitions(op, source, 1, collect(actual), parity));
    auto expected = source;
    expected[0] = 0;
    expected[69] = 1;
    EXPECT_EQ(actual, (std::vector<Transition>{{0, expected, {.25, .5}}}));
}

TEST(LocalTransitions, ValidatesEveryInputBeforeCallbacks) {
    Operator const valid{{2}, {{1., {{0, {0., 1., 1., 0.}}}, {}}}};
    std::vector<std::int64_t> const sources{0, 1};
    std::size_t callbacks = 0;
    auto visitor = [&](auto, auto const&, auto) { ++callbacks; };
    EXPECT_THROW(alps::ed::visit_local_transitions(valid, sources, 3, visitor), std::invalid_argument);
    EXPECT_THROW(alps::ed::visit_local_transitions(valid, sources,
        std::numeric_limits<std::size_t>::max(), visitor), std::invalid_argument);
    for (auto const& bad_source : {std::vector<std::int64_t>{-1, 0}, {0, 2}})
        EXPECT_THROW(alps::ed::visit_local_transitions(valid, bad_source, 2, visitor), std::invalid_argument);
    for (auto const& bad_parity : {std::vector<std::vector<double>>{{1.}}, {{1., -1.}, {1., -1.}}})
        EXPECT_THROW(alps::ed::visit_local_transitions(valid, sources, 2, visitor, bad_parity),
                     std::invalid_argument);
    for (int fault = 0; fault < 5; ++fault) {
        auto bad = valid;
        bad.terms.push_back(valid.terms.front());
        if (fault == 0) bad.site_dims[0] = 0;
        if (fault == 1) bad.terms.back().factors[0].site = 1;
        if (fault == 2) bad.terms.back().factors[0].matrix = {1.};
        if (fault == 3) bad.terms.back().factors.push_back(bad.terms.back().factors.front());
        if (fault == 4) bad.terms.back().jordan_wigner = {1};
        EXPECT_THROW(alps::ed::visit_local_transitions(bad, sources, 2, visitor), std::invalid_argument);
    }
    EXPECT_EQ(callbacks, 0);
}

TEST(LocalTransitions, PropagatesCallbackExceptionsAndRejectsEmptyVisitor) {
    Operator const op{{2}, {{1., {}, {}}, {2., {}, {}}}};
    std::vector<std::int64_t> const sources{0, 1};
    struct Failure {};
    std::size_t callbacks = 0;
    auto visitor = [&](auto, auto const&, auto) { ++callbacks; throw Failure{}; };
    EXPECT_THROW(alps::ed::visit_local_transitions(op, sources, 2, visitor), Failure);
    EXPECT_EQ(callbacks, 1);
    EXPECT_THROW(alps::ed::visit_local_transitions(op, sources, 2, {}), std::invalid_argument);
    EXPECT_THROW(alps::ed::visit_local_transitions(op, {}, 0, {}), std::invalid_argument);
}

TEST(LocalTransitions, HandlesVacuumEmptyBatchAndZeroParity) {
    Operator op{{}, {{{0., 2.}, {}, {}}}};
    std::vector<Transition> actual;
    EXPECT_TRUE(alps::ed::visit_local_transitions(op, {}, 1, collect(actual)));
    EXPECT_EQ(actual, (std::vector<Transition>{{0, {}, {0., 2.}}}));
    actual.clear();
    EXPECT_FALSE(alps::ed::visit_local_transitions(op, {}, 0, collect(actual)));
    EXPECT_TRUE(actual.empty());
    op = {{1}, {{{0., 2.}, {}, {0}}}};
    std::vector<std::int64_t> const sources{0};
    EXPECT_FALSE(alps::ed::visit_local_transitions(op, sources, 1, collect(actual), {{0.}}));
    EXPECT_TRUE(actual.empty());
}
