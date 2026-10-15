// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// ALPS's own test of params::supply() and params::supplied_text(), local
// additions to ALPSCore's params; see the README.

#include <alps/params.hpp>
#include <alps/testing/unique_file.hpp>
#include <gtest/gtest.h>

#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// The error message has_missing() prints, or "" if there is none
std::string missing_message(alps::params const& p) {
    std::ostringstream out;
    p.has_missing(out);
    return out.str();
}

} // namespace

TEST(ParamsSupply, TextIsParsedWithTheDefinedType) {
    alps::params p;
    p.supply("N", "24").supply("BETA", "12.3456789").supply("NAME", "a b");
    p.define<int>("N", "a size")
     .define<double>("BETA", "inverse temperature")
     .define<std::string>("NAME", "a name");
    ASSERT_TRUE(p.ok()) << missing_message(p);
    EXPECT_EQ(24, p["N"].as<int>());
    EXPECT_EQ(12.3456789, p["BETA"].as<double>());
    EXPECT_EQ("a b", p["NAME"].as<std::string>());
    EXPECT_TRUE(p.supplied("N"));
    EXPECT_FALSE(p.defaulted("N"));
}

TEST(ParamsSupply, SuppliedTextOverridesTheDefault) {
    alps::params p;
    p.supply("N", "3");
    p.define<int>("N", 7, "a size").define<int>("M", 7, "another size");
    ASSERT_TRUE(p.ok()) << missing_message(p);
    EXPECT_EQ(3, p["N"].as<int>());
    EXPECT_EQ(7, p["M"].as<int>());
    EXPECT_TRUE(p.defaulted("M"));
}

// A non-integral or malformed number for an integer parameter is an error
// naming the parameter.
TEST(ParamsSupply, NonIntegralTextForAnIntegerIsAnError) {
    for (std::string const text : {"2.7", "4.0", "1e6", "12abc", ""}) {
        alps::params p;
        p.supply("N_TAU", text);
        p.define<int>("N_TAU", 1000, "time slices");
        EXPECT_TRUE(p.has_missing()) << "'" << text << "'";
        EXPECT_NE(std::string::npos, missing_message(p).find("N_TAU")) << "'" << text << "'";
    }
}

TEST(ParamsSupply, ExponentTextIsADouble) {
    alps::params p;
    p.supply("SWEEPS", "1e6");
    p.define<double>("SWEEPS", "sweeps");
    ASSERT_TRUE(p.ok()) << missing_message(p);
    EXPECT_EQ(1e6, p["SWEEPS"].as<double>());
}

TEST(ParamsSupply, BooleanSpellings) {
    std::map<std::string, bool> const spellings = {
        {"true", true}, {"1", true}, {"yes", true}, {"on", true},
        {"false", false}, {"0", false}, {"no", false}, {"off", false}};
    for (auto const& spelling : spellings) {
        alps::params p;
        p.supply("FLAG", spelling.first);
        p.define<bool>("FLAG", !spelling.second, "a flag");
        ASSERT_TRUE(p.ok()) << spelling.first << ": " << missing_message(p);
        EXPECT_EQ(spelling.second, p["FLAG"].as<bool>()) << spelling.first;
    }
}

TEST(ParamsSupply, VectorsAreCommaSeparated) {
    alps::params p;
    p.supply("L", "1,2,3").supply("X", "0.5,-1e-3");
    p.define<std::vector<int> >("L", "lengths").define<std::vector<double> >("X", "positions");
    ASSERT_TRUE(p.ok()) << missing_message(p);
    EXPECT_EQ(std::vector<int>({1, 2, 3}), p["L"].as<std::vector<int> >());
    EXPECT_EQ(std::vector<double>({0.5, -1e-3}), p["X"].as<std::vector<double> >());
}

TEST(ParamsSupply, DefinedParametersCannotBeSupplied) {
    alps::params p;
    p.define<int>("N", 1, "a size");
    EXPECT_THROW(p.supply("N", "2"), std::logic_error);
    EXPECT_EQ(1, p["N"].as<int>());
}

// Only supplied text is listed, never defaults or assigned values, and text
// stays listed after its parameter is defined.
TEST(ParamsSupply, SuppliedTextListsWhatWasSupplied) {
    alps::params p;
    p.supply("U", "2").supply("H", "0.1");
    p.define<double>("U", "interaction").define<int>("SITES", 1, "sites");
    p["ASSIGNED"] = 5;
    std::map<std::string, std::string> const expected = {{"H", "0.1"}, {"U", "2"}};
    EXPECT_EQ(expected, p.supplied_text());
    std::ostringstream unused;
    EXPECT_TRUE(p.has_unused(unused));
    EXPECT_NE(std::string::npos, unused.str().find("H"));
}

TEST(ParamsSupply, SuppliedTextSurvivesAnArchive) {
    alps::testing::unique_file file("params_supply.h5.", alps::testing::unique_file::REMOVE_AFTER);
    alps::params p;
    p.supply("N", "24").supply("H", "0.1");
    p.define<int>("N", "a size");
    {
        alps::hdf5::archive ar(file.name(), "w");
        ar["/parameters"] << p;
    }
    alps::params restored;
    {
        alps::hdf5::archive ar(file.name(), "r");
        ar["/parameters"] >> restored;
    }
    EXPECT_EQ(p.supplied_text(), restored.supplied_text());
    EXPECT_EQ(24, restored["N"].as<int>());
}
