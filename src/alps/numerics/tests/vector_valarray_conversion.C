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

#include <alps/numeric/vector_valarray_conversion.hpp>
#include <gtest/gtest.h>
#include <valarray>
#include <vector>

TEST(VectorValarrayConversion, ValarrayToVectorPreservesSizeAndOrder)
{
    std::valarray<double> source(10);
    for (std::size_t i = 0; i < source.size(); ++i) source[i] = i;
    const auto result = alps::numeric::valarray2vector<double>(source);
    ASSERT_EQ(result.size(), source.size());
    for (std::size_t i = 0; i < source.size(); ++i)
        EXPECT_EQ(result[i], double(i)) << "index " << i;
}

TEST(VectorValarrayConversion, VectorToValarrayPreservesSizeAndOrder)
{
    std::vector<double> source;
    for (int i = 0; i < 10; ++i) source.push_back(10 - i);
    const auto result = alps::numeric::vector2valarray<double>(source);
    ASSERT_EQ(result.size(), source.size());
    for (std::size_t i = 0; i < source.size(); ++i)
        EXPECT_EQ(result[i], 10. - i) << "index " << i;
}

TEST(VectorValarrayConversion, EmptyContainersRemainEmpty)
{
    EXPECT_TRUE(alps::numeric::valarray2vector<double>(std::valarray<double>()).empty());
    EXPECT_EQ(alps::numeric::vector2valarray<double>(std::vector<double>()).size(), 0u);
    const auto converted = alps::numeric::vector2valarray<int, double>(std::vector<int>());
    EXPECT_EQ(converted.size(), 0u);
}

TEST(VectorValarrayConversion, ConvertingElementTypePreservesSizeAndOrder)
{
    const std::vector<int> source{3, -5, 8};
    const auto converted = alps::numeric::vector2valarray<int, double>(source);
    ASSERT_EQ(converted.size(), source.size());
    for (std::size_t i = 0; i < source.size(); ++i)
        EXPECT_DOUBLE_EQ(converted[i], static_cast<double>(source[i]));
}
