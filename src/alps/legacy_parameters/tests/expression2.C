/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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
#include <sstream>
#include <stdexcept>

TEST(ExpressionSimplification, CombinesLikeTermsAndPreservesValues) {
  const char *input[] = {"3*(a*b)*2+5*(x*y)^2*x*3", "3*a*b+5*a-3*a*b+2.5*a",
                         "3*a*b+5*a+3*a*b+2.5*a"};
  const char *canonical[] = {"15 * (x * y)^2 * x + 6 * a * b", "7.5 * a", "7.5 * a + 6 * a * b"};
  alps::Parameters params;
  params["a"] = 2;
  params["b"] = 3;
  params["x"] = 4;
  params["y"] = 5;
  for (int i = 0; i < 3; ++i) {
    SCOPED_TRACE(input[i]);
    alps::Expression expression(input[i]);
    const auto before = alps::evaluate<double>(expression, alps::ParameterEvaluator(params));
    expression.simplify();
    EXPECT_DOUBLE_EQ(alps::evaluate<double>(expression, alps::ParameterEvaluator(params)), before);
    std::ostringstream rendered;
    rendered << expression;
    EXPECT_EQ(rendered.str(), canonical[i]);
  }
}

// Simplifying used to turn x/(2*y) into 2*x*y, and printing 1/x gave x. The
// model code prints simplified terms and parses them again, so both the value
// and the printed text must keep the division.
TEST(ExpressionSimplification, KeepsDivisionBySymbols) {
  alps::Parameters params;
  params["x"] = 3;
  params["y"] = 5;
  for (const char *input : {"1/x", "y/x", "x/(2*y)", "x*x/(2*y)", "-1/x", "y/2"}) {
    SCOPED_TRACE(input);
    alps::Expression expression(input);
    const auto expected = alps::evaluate<double>(expression, alps::ParameterEvaluator(params));
    expression.simplify();
    EXPECT_DOUBLE_EQ(alps::evaluate<double>(expression, alps::ParameterEvaluator(params)), expected);
    std::ostringstream rendered;
    rendered << expression;
    alps::Expression reparsed(rendered.str());
    EXPECT_DOUBLE_EQ(alps::evaluate<double>(reparsed, alps::ParameterEvaluator(params)), expected)
        << rendered.str();
  }
}

// Evaluating a function no evaluator knows, such as the operator Sz(j) left in
// a coefficient, used to recurse until the stack overflowed (issue #385).
TEST(ExpressionEvaluation, RefusesUnknownFunction) {
  for (const char *input : {"Sz(j)", "2*Sz(j)", "exchange(i,j)"}) {
    SCOPED_TRACE(input);
    alps::expression::Expression<double> expression(input);
    EXPECT_FALSE(expression.can_evaluate());
    EXPECT_THROW(expression.value(), std::runtime_error);
  }
}

// The guard above must still accept every two-argument function the evaluator
// computes: atan2 with random evaluation off, and normal_random(mean, sigma).
TEST(ExpressionEvaluation, EvaluatesTwoArgumentFunctions) {
  const alps::expression::Evaluator<double> no_random(false);
  alps::expression::Expression<double> angle("atan2(1,1)");
  EXPECT_TRUE(angle.can_evaluate(no_random));
  EXPECT_DOUBLE_EQ(angle.value(no_random), std::atan2(1., 1.));

  alps::expression::Expression<double> normal("normal_random(1,0)");
  EXPECT_TRUE(normal.can_evaluate());
  EXPECT_DOUBLE_EQ(normal.value(), 1.);
}
