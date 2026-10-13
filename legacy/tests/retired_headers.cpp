// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// Headers with no remaining users live only in the legacy layer until 4.0.
// They must keep compiling while they are installed. <alps/utility/factory.hpp>
// is the exception: it has not compiled since before 3.0.0 (a syntax error in
// factory::register_type), so no consumer can depend on it.
#include <gtest/gtest.h>

#include <alps/utility/generator.hpp>

TEST(LegacyRetiredHeaders, StillCompile)
{
    SUCCEED();
}
