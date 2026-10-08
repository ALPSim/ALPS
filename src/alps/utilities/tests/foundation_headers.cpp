// SPDX-License-Identifier: MIT

// This consumer must compile with the utilities/container include roots only.
// It exercises array storage and formatting without the numerical umbrella.
#include <gtest/gtest.h>
#include <alps/multi_array/multi_array.hpp>
#include <alps/ngs/short_print.hpp>
#include <alps/utility/resize.hpp>

#include <iostream>
#include <sstream>

TEST(FoundationHeaders, IndependentArrayStorageAndFormatting)
{
    alps::multi_array<int, 2> values(2, 3);
    for (std::size_t i = 0; i < values.num_elements(); ++i)
        values.data()[i] = static_cast<int>(i + 1);
    alps::multi_array<int, 2> copy;
    copy = values;
    ASSERT_EQ(copy.shape()[0], 2u);
    ASSERT_EQ(copy.shape()[1], 3u);
    EXPECT_EQ(copy[1][2], 6);
    copy[0][0] = 9;
    EXPECT_EQ(values[0][0], 1);
    alps::multi_array<double, 2> resized(1, 1);
    alps::resize_same_as(resized, values);
    EXPECT_EQ(resized.shape()[0], 2u);
    EXPECT_EQ(resized.shape()[1], 3u);
    std::ostringstream output;
    output << alps::short_print(values);
    EXPECT_EQ(output.str(), "[1,..6..,6]");
}
