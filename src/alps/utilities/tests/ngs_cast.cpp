// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/utilities/cast.hpp>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <string>

TEST(NgsCastFromString, ReadsWholeNumbers) {
    EXPECT_EQ(alps::cast<int>(std::string("42")), 42);
    EXPECT_EQ(alps::cast<int>(std::string("  -7 ")), -7);
    EXPECT_EQ(alps::cast<long>(std::string("2147483648")), 2147483648L);
    EXPECT_EQ(alps::cast<unsigned>(std::string("4294967295")), 4294967295u);
    EXPECT_EQ(alps::cast<double>(std::string("0.707106781186547")), 0.707106781186547);
    EXPECT_EQ(alps::cast<double>(std::string("-2.5E+2")), -250.);
    EXPECT_EQ(alps::cast<double>(std::string(" 1e-3\n")), 1e-3);
    EXPECT_EQ(alps::cast<float>(std::string("0.5")), 0.5f);
    EXPECT_EQ(alps::cast<long double>(std::string("0.25")), 0.25L);
    EXPECT_EQ(alps::cast<double>(std::string("inf")), std::numeric_limits<double>::infinity());
}

// Decimal text for an integer keeps its historical integer value.
TEST(NgsCastFromString, IntegerFromDecimalText) {
    EXPECT_EQ(alps::cast<int>(std::string("4.0")), 4);
    EXPECT_EQ(alps::cast<int>(std::string("4.5")), 4);
    EXPECT_EQ(alps::cast<int>(std::string("-4.5")), -4);
    EXPECT_EQ(alps::cast<int>(std::string("1e3")), 1000);
    EXPECT_EQ(alps::cast<unsigned>(std::string("1e3")), 1000u);
}

TEST(NgsCastFromString, RejectsTextThatIsNotANumber) {
    for (char const * text : {"", "  ", "6x", "2*L", "U/2", "Hirsch-Fye", "1,5", "1 2", "nan?"}) {
        EXPECT_THROW(alps::cast<int>(std::string(text)), std::invalid_argument) << text;
        EXPECT_THROW(alps::cast<unsigned>(std::string(text)), std::invalid_argument) << text;
        EXPECT_THROW(alps::cast<double>(std::string(text)), std::invalid_argument) << text;
        EXPECT_THROW(alps::cast<float>(std::string(text)), std::invalid_argument) << text;
    }
}

TEST(NgsCastFromString, RejectsIntegersOutOfRange) {
    EXPECT_THROW(alps::cast<int>(std::string("2147483648")), std::out_of_range);
    EXPECT_THROW(alps::cast<int>(std::string("1e10")), std::out_of_range);
    EXPECT_THROW(alps::cast<unsigned>(std::string("-1")), std::out_of_range);
    EXPECT_THROW(alps::cast<short>(std::string("40000")), std::out_of_range);
}

TEST(NgsCastFromString, ReadsBooleanWordsAndNumbers) {
    for (char const * text : {"true", "True", " TRUE ", "1", "2"})
        EXPECT_TRUE(alps::cast<bool>(std::string(text))) << text;
    for (char const * text : {"false", "False", "0", "0.0"})
        EXPECT_FALSE(alps::cast<bool>(std::string(text))) << text;
    for (char const * text : {"", "yes", "truthy"})
        EXPECT_THROW(alps::cast<bool>(std::string(text)), std::invalid_argument) << text;
}

TEST(NgsCastFromString, RoundTripsThroughText) {
    EXPECT_EQ(alps::cast<double>(alps::cast<std::string>(0.707106781186547)), 0.707106781186547);
    EXPECT_EQ(alps::cast<int>(alps::cast<std::string>(-123)), -123);
    EXPECT_EQ(alps::cast<bool>(alps::cast<std::string>(true)), true);
}
