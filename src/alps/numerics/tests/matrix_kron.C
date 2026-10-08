/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2026 by the ALPS collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Layout test for alps::numeric::kron().
//
// Standard Kronecker product layout:
//   (A (x) B)(i1*rows(B) + i2, j1*cols(B) + j2) = A(i1,j1) * B(i2,j2)
// i.e. the INNER factor B's dimensions stride the output. The historical
// implementation used A's dimensions in the index formula, which coincides
// with the standard layout only for equal-size operands; for unequal shapes
// it scattered entries (and indexed out of bounds when A is the larger
// factor). The equal-size arm below doubles as a no-behavior-change check
// for callers passing equal-size square operands.

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>

#include <gtest/gtest.h>
#include <array>
#include <cstddef>

class KroneckerLayout : public ::testing::TestWithParam<std::array<std::size_t, 4>> {};

TEST_P(KroneckerLayout, InnerFactorStridesOutput)
{
    const auto [ar, ac, br, bc] = GetParam();
    alps::numeric::matrix<double> a(ar, ac), b(br, bc);
    for (std::size_t i = 0; i < ar; ++i)
        for (std::size_t j = 0; j < ac; ++j) a(i, j) = 1.25 + 10. * i + j;
    for (std::size_t i = 0; i < br; ++i)
        for (std::size_t j = 0; j < bc; ++j) b(i, j) = -3.5 + 10. * i + j;

    const auto result = alps::numeric::kron(a, b);
    ASSERT_EQ(num_rows(result), ar * br);
    ASSERT_EQ(num_cols(result), ac * bc);
    for (std::size_t i1 = 0; i1 < ar; ++i1)
        for (std::size_t j1 = 0; j1 < ac; ++j1)
            for (std::size_t i2 = 0; i2 < br; ++i2)
                for (std::size_t j2 = 0; j2 < bc; ++j2)
                    EXPECT_EQ(result(i1 * br + i2, j1 * bc + j2), a(i1, j1) * b(i2, j2))
                        << "output row " << i1 * br + i2 << ", column " << j1 * bc + j2;
}

INSTANTIATE_TEST_SUITE_P(RectangularAndSquare, KroneckerLayout, ::testing::Values(
    std::array<std::size_t, 4>{3, 3, 3, 3},
    std::array<std::size_t, 4>{2, 3, 4, 2},
    std::array<std::size_t, 4>{4, 2, 2, 5},
    std::array<std::size_t, 4>{1, 1, 3, 2}));
