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

#include <alps/ngs/params.hpp>

#include <gtest/gtest.h>
#include <vector>

TEST(ParamsAssignment, NativeTypesAndNonInsertingLookup) {

    alps::params parms;
    parms["char"] = static_cast<char>(1);
    parms["signed char"] = static_cast<signed char>(1);
    parms["unsigned char"] = static_cast<unsigned char>(1);
    parms["short"] = static_cast<short>(1);
    parms["unsigned short"] = static_cast<unsigned short>(1);
    parms["int"] = static_cast<int>(1);
    parms["unsigned"] = static_cast<unsigned>(1);
    parms["long"] = static_cast<long>(1);
    parms["unsigned long"] = static_cast<unsigned long>(1);
    parms["long long"] = static_cast<long long>(1);
    parms["unsigned long long"] = static_cast<unsigned long long>(1);
    parms["float"] = static_cast<float>(1);
    parms["double"] = static_cast<double>(1);
    parms["long double"] = static_cast<long double>(1);
    parms["bool"] = static_cast<bool>(1);
    std::vector<bool> const bool_vector{true, false, true};
    parms["std::vector<bool>"] = bool_vector;
    parms["std::string"] = std::string("asdf");

    EXPECT_EQ(parms["std::vector<bool>"].cast<std::vector<bool> >(), bool_vector);
    ASSERT_NE(parms.find("int"), nullptr);
    EXPECT_EQ(parms.find("int")->cast<int>(), 1);
    EXPECT_EQ(parms.find("missing"), nullptr);

    for (const auto& name : {"char", "signed char", "unsigned char", "short", "unsigned short",
                             "int", "unsigned", "long", "unsigned long", "long long",
                             "unsigned long long", "float", "double", "long double", "bool"}) {
        EXPECT_EQ(parms[name].cast<int>(), 1) << name;
    }
    EXPECT_EQ(parms["std::string"].cast<std::string>(), "asdf");
}
