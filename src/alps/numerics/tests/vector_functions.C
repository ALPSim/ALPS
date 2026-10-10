/****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id: nobinning.h 3520 2009-12-11 16:49:53Z gamperl $ */

#include <alps/numeric/vector_functions.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <string>
#include <vector>

namespace {
void expect_elements(const std::vector<double>& actual, double expected,
                     const char* operation)
{
    SCOPED_TRACE(operation);
    ASSERT_EQ(actual.size(), 10u);
    for (std::size_t i = 0; i < actual.size(); ++i)
        // These element-wise formulas perform only a handful of double operations.
        EXPECT_NEAR(actual[i], expected, 2e-14) << "index " << i;
}
}

TEST(VectorFunctions, ArithmeticAndTranscendentalsMatchScalarReferences)
{
    using namespace alps::numeric;
    std::vector<double> vecA(10, 0.81), vecB(10, 1.21);
    expect_elements(+vecA, +0.81, "+vecA");
    expect_elements(-vecA, -0.81, "-vecA");
    expect_elements(vecA + vecB, 0.81 + 1.21, "vecA + vecB");
    expect_elements(vecA - vecB, 0.81 - 1.21, "vecA - vecB");
    expect_elements(vecA * vecB, 0.81 * 1.21, "vecA * vecB");
    expect_elements(vecA / vecB, 0.81 / 1.21, "vecA / vecB");
    expect_elements(vecA + 1., 0.81 + 1., "vecA + 1.");
    expect_elements(vecA - 1., 0.81 - 1., "vecA - 1.");
    expect_elements(vecA * 1., 0.81 * 1., "vecA * 1.");
    expect_elements(vecA / 1., 0.81 / 1., "vecA / 1.");
    expect_elements(1. + vecA, 1. + 0.81, "1. + vecA");
    expect_elements(1. - vecA, 1. - 0.81, "1. - vecA");
    expect_elements(1. * vecA, 1. * 0.81, "1. * vecA");
    expect_elements(1. / vecA, 1. / 0.81, "1. / vecA");
    expect_elements(pow(vecA,2.71), std::pow(0.81,2.71), "pow(vecA,2.71)");
    expect_elements(sq(vecA), (0.81 * 0.81), "sq(vecA)");
    expect_elements(sqrt(vecA), std::sqrt(0.81), "sqrt(vecA)");
    expect_elements(cb(vecA), (0.81 * 0.81 * 0.81), "cb(vecA)");
    expect_elements(cbrt(vecA), std::cbrt(0.81), "cbrt(vecA)");
    expect_elements(exp(vecA), std::exp(0.81), "exp(vecA)");
    expect_elements(log(vecA), std::log(0.81), "log(vecA)");
    expect_elements(sin(vecA), std::sin(0.81), "sin(vecA)");
    expect_elements(cos(vecA), std::cos(0.81), "cos(vecA)");
    expect_elements(tan(vecA), std::tan(0.81), "tan(vecA)");
    expect_elements(asin(vecA), std::asin(0.81), "asin(vecA)");
    expect_elements(acos(vecA), std::acos(0.81), "acos(vecA)");
    expect_elements(atan(vecA), std::atan(0.81), "atan(vecA)");
    expect_elements(sinh(vecA), std::sinh(0.81), "sinh(vecA)");
    expect_elements(cosh(vecA), std::cosh(0.81), "cosh(vecA)");
    expect_elements(tanh(vecA), std::tanh(0.81), "tanh(vecA)");
    expect_elements(asinh(vecA), std::asinh(0.81), "asinh(vecA)");
    expect_elements(acosh(vecB), std::acosh(1.21), "acosh(vecB)");
    expect_elements(atanh(vecA), std::atanh(0.81), "atanh(vecA)");
    expect_elements(sqrt(vecA*vecA + vecB*vecB), std::sqrt(0.81*0.81 + 1.21*1.21), "sqrt(vecA*vecA + vecB*vecB)");
}
