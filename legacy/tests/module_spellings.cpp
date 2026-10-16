// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// Old alps/utility spellings of headers owned by ALPS::cli and ALPS::numerics
// forward to the owner's prefix.
#include <gtest/gtest.h>

#include <alps/utility/cli.hpp>
#include <alps/utility/numeric_cast.hpp>

#include <complex>

TEST(LegacyModuleSpellings, UtilityHeadersForwardToOwners)
{
    EXPECT_TRUE(alps::cli_is_master());
    EXPECT_EQ(alps::numeric_cast<double>(std::complex<double>(2.5, 0.0)), 2.5);
}
