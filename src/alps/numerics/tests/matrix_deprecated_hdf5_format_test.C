/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2013        by Michele Dolfi <dolfim@phys.ethz.ch>,               *
 *                              Andreas Hehn <hehn@phys.ethz.ch>                   *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <alps/hdf5/matrix.hpp>
#include <alps/numeric/matrix.hpp>
#include <alps/hdf5.hpp>
#include <gtest/gtest.h>
#include <filesystem>
#include <cmath>

TEST(MatrixCompatibility, ReadsHistoricalHdf5Layout)
{
    const auto filename = std::filesystem::path(ALPS_TEST_DATA_DIR) /
        "matrix_deprecated_hdf5_format_test.h5";
    ASSERT_TRUE(std::filesystem::exists(filename)) << filename;
    alps::numeric::matrix<double> matrix;
    alps::hdf5::archive archive(filename.string(), "r");
    archive["/matrix_old_hdf5_format"] >> matrix;
    ASSERT_EQ(num_rows(matrix), 20u);
    ASSERT_EQ(num_cols(matrix), 3u);
    EXPECT_EQ(matrix.capacity().first, 30u);
    EXPECT_EQ(matrix.capacity().second, 3u);
    // Historical stdout fixture recorded six significant digits. Preserve those
    // independent reference values with a bound matching that rounding precision.
    const double expected[20][3] = {
        {-0.548267, -0.5067, 0.0886095},
        {-0.353845, 0.0520657, -0.66079},
        {0.0132811, 0.108343, -0.22786},
        {0.114659, -0.347019, 0.264671},
        {0.0795533, -0.0638698, -0.26175},
        {0.00714078, -0.0156018, -0.057123},
        {0.00274422, -0.0883628, -0.0458315},
        {-0.00984287, -0.028685, 0.0536734},
        {0.000211855, 0.0562105, 0.0482957},
        {-0.00180964, -0.00111478, -0.0424608},
        {-0.000332624, -0.00743965, 0.00493614},
        {0.662257, 0.0202139, -0.131487},
        {-0.193568, 0.608751, -0.0471581},
        {0.070073, -0.165623, -0.457068},
        {0.0038202, 0.132761, -0.0140661},
        {-0.00172488, -0.0172448, -0.147131},
        {0.000666784, 0.0220719, -0.0260034},
        {0.000111599, -0.0207896, -0.072886},
        {0.270933, -0.370897, -0.298945},
        {-0.000490229, -0.193364, -0.0836068},
    };
    for (std::size_t row = 0; row < 20; ++row)
        for (std::size_t column = 0; column < 3; ++column) {
            const double reference = expected[row][column];
            const double tolerance = 0.500001 * std::pow(10., std::floor(std::log10(std::abs(reference))) - 5);
            EXPECT_NEAR(matrix(row, column), reference, tolerance)
                << "row " << row << ", column " << column;
        }
}
