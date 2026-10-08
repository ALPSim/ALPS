/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2013 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// test for merging two empty histogram observables

#include <alps/alea.h>
#include <gtest/gtest.h>

TEST(Histogram, MergingEmptyObservablesPreservesRangeAndZeroCounts)
{
    alps::ObservableSet first, second;
    first << alps::IntHistogramObservable("histogram", 0, 10);
    second << alps::IntHistogramObservable("histogram", 0, 10);
    first << second;
    const alps::IntHistogramObsevaluator result(first["histogram"]);
    EXPECT_EQ(result.count(), 0u);
    ASSERT_EQ(result.size(), 10u);
    EXPECT_EQ(result.min(), 0);
    EXPECT_EQ(result.max(), 10);
    EXPECT_EQ(result.stepsize(), 1);
    for (std::size_t i = 0; i < result.size(); ++i) EXPECT_EQ(result[i], 0u);
}
