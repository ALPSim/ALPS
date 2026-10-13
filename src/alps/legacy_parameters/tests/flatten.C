/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2002 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

TEST(ExpressionFlattening, DistributesProductsAndPreservesValues) {
  const char *input[] = {"(a+b)*(c+d)", "(-a+b)/(c+d)", "sin(5*(x+y))", "(a+b)^2*(c+d)*(e+f)"};
  const char *canonical[] = {"(a) * (c) + (a) * (d) + (b) * (c) + (b) * (d)",
                             "( - a) / (c + d) + (b) / (c + d)", "sin(5 * (x) + 5 * (y))",
                             "(a + b)^2 * (c) * (e) + (a + b)^2 * (c) * (f) + (a + b)^2 * (d) * "
                             "(e) + (a + b)^2 * (d) * (f)"};
  alps::Parameters params;
  for (const auto *key : {"a", "b", "c", "d", "e", "f", "x", "y"})
    params[key] = 2;
  for (int i = 0; i < 4; ++i) {
    SCOPED_TRACE(input[i]);
    alps::Expression expression(input[i]);
    const auto before = alps::evaluate<double>(expression, alps::ParameterEvaluator(params));
    expression.flatten();
    EXPECT_NEAR(alps::evaluate<double>(expression, alps::ParameterEvaluator(params)), before,
                1e-12);
    std::ostringstream rendered;
    rendered << expression;
    EXPECT_EQ(rendered.str(), canonical[i]);
  }
}
