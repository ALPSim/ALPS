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
#include <vector>

TEST(ParamsOrdering, CheckpointPreservesInsertionOrderAndValues) {
    alps::testing::TemporaryDirectory temporary;
    const auto filename = (temporary.path() / "ordering.h5").string();
    alps::params expected;
    expected["a"] = 6;
    expected["x"] = 2;
    expected["b"] = 3;
    expected["w"] = 1;
    const std::vector<std::string> keys{"a", "x", "b", "w"};
    std::vector<std::string> before;
    for (const auto& entry : expected) before.push_back(entry.first);
    EXPECT_EQ(before, keys);
    { alps::hdf5::archive archive(filename, "w"); archive["/parameters"] << expected; }
    alps::params actual;
    { alps::hdf5::archive archive(filename, "r"); archive["/parameters"] >> actual; }
    std::vector<std::string> after;
    for (const auto& entry : actual) after.push_back(entry.first);
    ASSERT_EQ(after, keys);
    for (const auto& key : keys) EXPECT_EQ(actual[key].cast<int>(), expected[key].cast<int>()) << key;
}
