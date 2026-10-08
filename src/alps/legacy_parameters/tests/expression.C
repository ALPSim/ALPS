/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2006 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/expression.h>
#include <gtest/gtest.h>
#include <cmath>
#include <vector>

namespace {
alps::Parameters parameters() {
  alps::Parameters result;
  result["L"] = 10;
  result["T"] = 0.1;
  result["beta"] = "1/T";
  result["error"] = "error";
  return result;
}
struct Evaluation {
  const char *expression;
  double expected;
};
void PrintTo(const Evaluation &sample, std::ostream *output) { *output << sample.expression; }
class ExpressionEvaluation : public ::testing::TestWithParam<Evaluation> {};
TEST_P(ExpressionEvaluation, EvaluatesParametersAndFunctions) {
  const auto sample = GetParam();
  SCOPED_TRACE(sample.expression);
  const auto params = parameters();
  alps::ParameterEvaluator evaluator(params);
  alps::Expression expression(sample.expression);
  ASSERT_TRUE(expression.can_evaluate(evaluator));
  // Analytic references; bound rounding error in elementary operations.
  EXPECT_NEAR(alps::evaluate<double>(expression, evaluator), sample.expected, 1e-11);
  EXPECT_NEAR(alps::evaluate<double>(sample.expression, params), sample.expected, 1e-11);
}
INSTANTIATE_TEST_SUITE_P(
    HistoricalInputs, ExpressionEvaluation,
    ::testing::Values(Evaluation{"sqrt(4)", 2}, Evaluation{"3+5", 8}, Evaluation{"L", 10},
                      Evaluation{"T", .1}, Evaluation{"1/L", .1}, Evaluation{"1/T", 10},
                      Evaluation{"2*Pi/L", 2 * std::acos(-1.) / 10}, Evaluation{"L+T", 10.1},
                      Evaluation{"L+10", 20}, Evaluation{"L-T", 9.9}, Evaluation{"1/T+10", 20},
                      Evaluation{"L/(1/T+10)", .5}, Evaluation{"L*L", 100}, Evaluation{"beta", 10},
                      Evaluation{"sqrt(L)", std::sqrt(10.)},
                      Evaluation{"sin(2*Pi/L)", std::sin(2 * std::acos(-1.) / 10)},
                      Evaluation{"cos(2*Pi/L)", std::cos(2 * std::acos(-1.) / 10)},
                      Evaluation{"L^2", 100}, Evaluation{"(L+1)^3", 1331},
                      Evaluation{"(L+1)^-1*5", 5. / 11}, Evaluation{"sqrt(3.33)", std::sqrt(3.33)},
                      Evaluation{"L*T", 1}, Evaluation{"2*Pi*L", 20 * std::acos(-1.)}));

TEST(Expression, RejectsUnresolvedSelfReference) {
  const auto params = parameters();
  alps::ParameterEvaluator evaluator(params);
  EXPECT_FALSE(alps::Expression("error").can_evaluate(evaluator));
  EXPECT_FALSE(alps::can_evaluate(std::string("error"), params));
}

TEST(Expression, PreservesDefaultRandomSequence) {
  // Rounded references from expression.output. Each random function has its
  // own default engine; exercise calls in the same order as the legacy test.
  const auto params = parameters();
  const std::vector<Evaluation> samples{{"random()", .814724},
                                        {"random()", .135477},
                                        {"random()", .905792},
                                        {"integer_random(50)", 41},
                                        {"integer_random(50)", 6},
                                        {"integer_random(50)", 48},
                                        {"gaussian_random()", 2.47621},
                                        {"gaussian_random()", -.52885},
                                        {"gaussian_random(100,20)", 113.183}};
  for (const auto &sample : samples) {
    SCOPED_TRACE(sample.expression);
    EXPECT_NEAR(alps::evaluate<double>(sample.expression, params), sample.expected,
                std::abs(sample.expected) > 100 ? 5e-4 : 5e-6);
  }
}
} // namespace
