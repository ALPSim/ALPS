// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// A result file whose bytes are fixed: fixtures/legacy_observables.h5 was
// written once by ALEA (fixtures/write_legacy_observables.cpp), together with
// every number a reader got back from it then. Whatever reads result files
// must return exactly those numbers.
#include "fixtures/legacy_observables.hpp"

#include <alps/alea.h>
#include <alps/hdf5.hpp>

#include <gtest/gtest.h>

#include <cstdlib>
#include <fstream>
#include <string>

namespace {

const std::string fixtures = std::string(ALPS_TEST_SOURCE_DIR) + "/fixtures/";

statistics_test::numbers expected(const std::string & group) {
    statistics_test::numbers values;
    std::ifstream in(fixtures + "legacy_observables.expected");
    std::string key, value;
    while (in >> key >> value)
        if (key.compare(0, group.size() + 1, group + "/") == 0)
            values[key.substr(group.size() + 1)] = std::strtod(value.c_str(), nullptr);
    return values;
}

void expect_equal(const statistics_test::numbers & want, const statistics_test::numbers & got) {
    ASSERT_FALSE(want.empty());
    ASSERT_EQ(want.size(), got.size());
    for (const auto & value : want) {
        const auto found = got.find(value.first);
        ASSERT_NE(found, got.end()) << value.first;
        if (value.second != value.second)
            EXPECT_NE(found->second, found->second) << value.first;
        else
            EXPECT_EQ(value.second, found->second) << value.first;
    }
}

alps::ObservableSet load(const std::string & group) {
    alps::ObservableSet set;
    alps::hdf5::archive archive(fixtures + "legacy_observables.h5", "r");
    archive[group] >> set;
    return set;
}

}

TEST(FrozenFiles, AleaReadsARun)
{
    const alps::ObservableSet set = load("/run");
    expect_equal(expected("run"), statistics_test::describe_run(set));
    EXPECT_EQ(set.get<alps::SignedObservable<alps::RealObservable> >("O").sign_name(), "Sign");
}

TEST(FrozenFiles, AleaReadsMergedRuns)
{
    expect_equal(expected("merged"), statistics_test::describe_merged(load("/merged")));
}

// ALEA writes merged runs that hold a signed observable but does not read them
// back: the group lacks an attribute its own reader requires. It fails loudly.
// The group is in the file for the reader that replaces this one.
TEST(FrozenFiles, AleaRefusesMergedSignedRuns)
{
    EXPECT_THROW(load("/merged_signed"), alps::hdf5::path_not_found);
}

// The frozen numbers are not arbitrary: the stream is known, so the file can
// be held to exact values as well.
TEST(FrozenFiles, TheFileHoldsTheKnownStream)
{
    const statistics_test::numbers run = expected("run");
    const std::vector<double> x = statistics_test::ar1_series(1001, 0.7, 20000);
    double sum = 0;
    for (double value : x)
        sum += value;
    EXPECT_EQ(run.at("x/count"), 20000);
    EXPECT_EQ(run.at("x/mean"), sum / 20000);
    EXPECT_EQ(run.at("s/mean"), sum / 20000);
    EXPECT_EQ(run.at("v/mean/0"), sum / 20000);
    EXPECT_EQ(run.at("h/count"), 20000);
    // Signed average 1.5, error 0.0127 at this length.
    EXPECT_NEAR(run.at("O/mean"), 1.5, 0.05);
}
