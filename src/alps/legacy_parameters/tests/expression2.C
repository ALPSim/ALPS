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
#include <sstream>

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
