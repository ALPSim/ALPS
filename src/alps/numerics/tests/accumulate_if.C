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

#include <alps/numeric/accumulate_if.hpp>
#include <gtest/gtest.h>
#include <numeric>
#include <vector>

TEST(AccumulateIf, PredicateSelectsElements)
{
    std::vector<double> values(30);
    for (std::size_t i = 0; i < values.size(); ++i) values[i] = 0.1 * i;
    const auto result = alps::numeric::accumulate_if(
        values.begin(), values.end(), 0., [](double value) { return value <= 2.; });
    // Arithmetic series 0 + 0.1 + ... + 2.0; allow accumulated rounding.
    EXPECT_NEAR(result, 21., 1e-13);
}

TEST(AccumulateIf, CustomOperationAccumulatesSquares)
{
    std::vector<double> values(30);
    for (std::size_t i = 0; i < values.size(); ++i) values[i] = 0.1 * i;
    const auto result = alps::numeric::accumulate_if(
        values.begin(), values.end(), 0.,
        [](double sum, double value) { return sum + value * value; },
        [](double value) { return value > 1.; });
    // sum(k^2, k=11..29) / 100 = (29*30*59 - 10*11*21) / 600.
    EXPECT_NEAR(result, 81.7, 1e-12);
}

TEST(AccumulateIf, EmptyAndRejectedRangesPreserveInitialValue)
{
    const std::vector<double> empty;
    EXPECT_EQ(alps::numeric::accumulate_if(empty.begin(), empty.end(), 7.,
        [](double) { return true; }), 7.);
    const std::vector<double> values{1., 2.};
    EXPECT_EQ(alps::numeric::accumulate_if(values.begin(), values.end(), 7.,
        [](double) { return false; }), 7.);
}
