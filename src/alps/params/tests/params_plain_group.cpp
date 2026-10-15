// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
//
// ALPS's own test of loading ALPSCore's params from an HDF5 group without
// ALPSCore's attributes, a local addition; see the README.

#include <alps/hdf5/complex.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <alps/params/text_values.hpp>
#include <alps/testing/unique_file.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <complex>
#include <cstdint>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string missing_message(alps::params const& p) {
    std::ostringstream out;
    p.has_missing(out);
    return out.str();
}

alps::params load(std::string const& file) {
    alps::hdf5::archive ar(file, "r");
    alps::params p;
    ar["/parameters"] >> p;
    return p;
}

class ParamsPlainGroup : public ::testing::Test {
  protected:
    alps::testing::unique_file file{"params_plain_group.h5.", alps::testing::unique_file::REMOVE_AFTER};
};

} // namespace

TEST(ParamsToText, Numbers) {
    using alps::params_ns::to_text;
    EXPECT_EQ("24", to_text(24.0));
    EXPECT_EQ("-3", to_text(-3.0));
    EXPECT_EQ("1000000", to_text(1e6));
    EXPECT_EQ("0.1", to_text(0.1));
    EXPECT_EQ("-0.5", to_text(-0.5));
    EXPECT_EQ("12.3456789", to_text(12.3456789));
    EXPECT_EQ("0.3333333333333333", to_text(1.0 / 3.0));
    EXPECT_EQ("1e+20", to_text(1e20));
    EXPECT_EQ("1", to_text(true));
    EXPECT_EQ("0", to_text(false));
    EXPECT_EQ("-7", to_text(-7));
    EXPECT_EQ("5000000000", to_text(5000000000LL));
    EXPECT_EQ("1,2,3", to_text(std::vector<int>({1, 2, 3})));
    EXPECT_EQ("0.5,24", to_text(std::vector<double>({0.5, 24.0})));
    EXPECT_EQ("1,0", to_text(std::vector<bool>({true, false})));
}

// Values as pyalps, legacy alps::Parameters and the NGS params wrote them
TEST_F(ParamsPlainGroup, SuppliesEachValueAsText) {
    {
        alps::hdf5::archive ar(file.name(), "w");
        ar["/parameters/N_TAU"] << 1000;
        ar["/parameters/BETA"] << 10.0;
        ar["/parameters/MU"] << 0.123456789;
        ar["/parameters/SWEEPS"] << 1e6;
        ar["/parameters/BIG"] << std::int64_t(5000000000LL);
        ar["/parameters/WIDE"] << 7u;
        ar["/parameters/FLAG"] << true;
        ar["/parameters/DELTA"] << std::string("delta.dat");
        ar["/parameters/L"] << std::vector<int>({4, 4});
        ar["/parameters/X"] << std::vector<double>({0.5, 1.0});
        ar["/parameters/NAMES"] << std::vector<std::string>({"a", "b"});
        // pyalps: a bool is an int8 tagged as such
        ar["/parameters/MEASURE"] << static_cast<signed char>(1);
        ar["/parameters/MEASURE/@__alps_type__"] << std::string("bool");
        // pyalps: a list of mixed elements is a group of numbered scalars
        ar["/parameters/MIXED/0"] << 1;
        ar["/parameters/MIXED/1"] << 2.5;
    }
    alps::params p = load(file.name());
    const std::map<std::string, std::string> expected = {
        {"N_TAU", "1000"}, {"BETA", "10"}, {"MU", "0.123456789"}, {"SWEEPS", "1000000"},
        {"BIG", "5000000000"}, {"WIDE", "7"}, {"FLAG", "1"}, {"DELTA", "delta.dat"},
        {"L", "4,4"}, {"X", "0.5,1"}, {"NAMES", "a,b"}, {"MEASURE", "1"}, {"MIXED", "1,2.5"}};
    EXPECT_EQ(expected, p.supplied_text());

    p.define<int>("N_TAU", "time slices")
     .define<double>("BETA", "inverse temperature")
     .define<int>("SWEEPS", "sweeps")
     .define<long>("BIG", "a large number")
     .define<bool>("FLAG", false, "a flag")
     .define<int>("MEASURE", "a flag read as an integer")
     .define<std::vector<int> >("L", "lengths")
     .define<std::vector<double> >("MIXED", "a mixed list");
    ASSERT_TRUE(p.ok()) << missing_message(p);
    EXPECT_EQ(1000, p["N_TAU"].as<int>());
    EXPECT_EQ(10.0, p["BETA"].as<double>());
    EXPECT_EQ(1000000, p["SWEEPS"].as<int>());
    EXPECT_EQ(5000000000L, p["BIG"].as<long>());
    EXPECT_TRUE(p["FLAG"].as<bool>());
    EXPECT_EQ(1, p["MEASURE"].as<int>());
    EXPECT_EQ(std::vector<int>({4, 4}), p["L"].as<std::vector<int> >());
    EXPECT_EQ(std::vector<double>({1.0, 2.5}), p["MIXED"].as<std::vector<double> >());
    EXPECT_FALSE(p.is_restored());
}

// A non-integral number for an integer parameter is an error naming it.
TEST_F(ParamsPlainGroup, NonIntegralNumberForAnInteger) {
    {
        alps::hdf5::archive ar(file.name(), "w");
        ar["/parameters/N_TAU"] << 2.7;
    }
    alps::params p = load(file.name());
    p.define<int>("N_TAU", 1000, "time slices");
    EXPECT_TRUE(p.has_missing());
    EXPECT_NE(std::string::npos, missing_message(p).find("N_TAU"));
}

TEST_F(ParamsPlainGroup, ComplexValuesAreAnErrorNamingThem) {
    {
        alps::hdf5::archive ar(file.name(), "w");
        ar["/parameters/Z"] << std::complex<double>(1.0, 2.0);
    }
    try {
        load(file.name());
        FAIL() << "a complex value loaded";
    } catch (std::runtime_error const& error) {
        EXPECT_NE(std::string::npos, std::string(error.what()).find("'Z'")) << error.what();
    }
}

// A group saved by ALPSCore's params keeps its own format.
TEST_F(ParamsPlainGroup, ALPSCoreGroupsLoadAsBefore) {
    alps::params saved;
    saved.supply("N", "3");
    saved.define<int>("N", "a size").define<double>("BETA", 2.0, "inverse temperature");
    {
        alps::hdf5::archive ar(file.name(), "w");
        ar["/parameters"] << saved;
    }
    alps::params p = load(file.name());
    EXPECT_TRUE(p == saved);
    EXPECT_TRUE(p.defined("BETA"));
    EXPECT_TRUE(p.defaulted("BETA"));
}

TEST(ParamsSuppliedValue, NumbersListsAndText) {
    using alps::params_ns::supplied_value;
    EXPECT_EQ(10L, supplied_value("N", "10").as<long>());
    EXPECT_EQ(-3L, supplied_value("N", "-3").as<long>());
    EXPECT_EQ(2.25, supplied_value("T", "2.25").as<double>());
    EXPECT_EQ(1e6, supplied_value("T", "1e6").as<double>());
    EXPECT_EQ(std::vector<long>({4, 4}), supplied_value("L", "4,4").as<std::vector<long> >());
    EXPECT_EQ(std::vector<double>({1.0, 2.5}), supplied_value("X", "1,2.5").as<std::vector<double> >());
    for (const char* text : {"", "chain lattice", "0x10", "1,a", "1,,2", "5 ", "true"}) {
        EXPECT_EQ(text, supplied_value("S", text).as<std::string>()) << text;
    }
}

// Supplied parameters the program never defines are datasets too, as ALPS's
// other writers left them, and load back as text.
TEST_F(ParamsPlainGroup, SuppliedOnlyValuesAreWrittenAndStayText) {
    alps::params saved;
    saved.supply("L", "10").supply("T", "2.25").supply("LATTICE", "chain").supply("X", "1,2.5")
         .supply("N", "3");
    saved.define<int>("N", "a size");
    {
        alps::hdf5::archive ar(file.name(), "w");
        ar["/parameters"] << saved;
    }
    {
        alps::hdf5::archive ar(file.name(), "r");
        std::vector<std::string> children = ar.list_children("/parameters");
        std::sort(children.begin(), children.end());
        EXPECT_EQ(std::vector<std::string>({"L", "LATTICE", "N", "T", "X"}), children);
        long length = 0;
        ar["/parameters/L"] >> length;
        EXPECT_EQ(10, length);
        EXPECT_TRUE(ar.is_datatype<long>("/parameters/L"));
        EXPECT_TRUE(ar.is_datatype<double>("/parameters/T"));
        std::vector<double> x;
        ar["/parameters/X"] >> x;
        EXPECT_EQ(std::vector<double>({1.0, 2.5}), x);
    }
    alps::params p = load(file.name());
    EXPECT_TRUE(p == saved);
    EXPECT_FALSE(p.exists("L"));
    // An integer-looking value parses as the type the program declares
    p.define<double>("L", "length").define<std::string>("T", "text");
    ASSERT_TRUE(p.ok()) << missing_message(p);
    EXPECT_EQ(10.0, p["L"].as<double>());
    EXPECT_EQ("2.25", p["T"].as<std::string>());
}

TEST_F(ParamsPlainGroup, EmptyParamsRoundTrip) {
    alps::params saved;
    {
        alps::hdf5::archive ar(file.name(), "w");
        ar["/parameters"] << saved;
    }
    EXPECT_TRUE(load(file.name()) == saved);
}
