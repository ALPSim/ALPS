/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                              Matthias Troyer <troyer@comp-phys.org>             *
 *                                                                                 *
 * ALPS Project: https://alps.comp-phys.org/                                       *
 * SPDX-License-Identifier: MIT                                                    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>
#include <alps/ngs/params.hpp>

TEST(ParamsDefault, MissingValueUsesFallbackWithoutInserting) {
    alps::params parameters;
    EXPECT_EQ(std::string(parameters["missing"] | "substitution_string"), "substitution_string");
    EXPECT_EQ(parameters["missing"] | 42, 42);
    EXPECT_FALSE(parameters.defined("missing"));
    parameters["present"] = 7;
    EXPECT_EQ(parameters["present"] | 42, 7);
}
