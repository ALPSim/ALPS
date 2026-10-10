/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2013 by Andreas Hehn <hehn@phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include <alps/utility/bitops.hpp>
#include <boost/cstdint.hpp>

#include <gtest/gtest.h>

TEST(BitOperations, CountsHighBitsInBothWidths) {
    const boost::uint32_t ui = boost::uint32_t(5) << (sizeof(boost::uint32_t)*8-4);
    const boost::uint64_t uli = boost::uint64_t(5) << (sizeof(boost::uint64_t)*8-4);
    EXPECT_EQ(alps::popcnt(ui), 2);
    EXPECT_EQ(alps::popcnt(uli), 2);
}
