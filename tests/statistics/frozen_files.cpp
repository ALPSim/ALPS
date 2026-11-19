// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// A result file whose bytes are fixed: fixtures/legacy_observables.h5 was
// written once by the ALEA of 3.0 (the writer went with it), together with
// every number its reader got back from it then. Whatever reads result files
// must return exactly those numbers.
#include "streams.hpp"

#include <alps/alea/observables_io.hpp>
#include <alps/hdf5.hpp>

#include <gtest/gtest.h>

#include <cstdlib>
#include <fstream>
#include <map>
#include <string>

namespace {

typedef std::map<std::string, double> numbers;

const std::string fixtures = std::string(ALPS_TEST_SOURCE_DIR) + "/fixtures/";

numbers expected(const std::string & group) {
    numbers values;
    std::ifstream in(fixtures + "legacy_observables.expected");
    std::string key, value;
    while (in >> key >> value)
        if (key.compare(0, group.size() + 1, group + "/") == 0)
            values[key.substr(group.size() + 1)] = std::strtod(value.c_str(), nullptr);
    return values;
}

void expect_equal(const numbers & want, const numbers & got) {
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

}

// The frozen numbers are not arbitrary: the stream is known, so the file can
// be held to exact values as well.
TEST(FrozenFiles, TheFileHoldsTheKnownStream)
{
    const numbers run = expected("run");
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

// The reader of the observable layer returns the stored numbers of the same
// file. It does not offer the binning analysis of 3.0 files (pyalps reads that
// itself), so the two per-level errors are not compared.
TEST(FrozenFiles, ObservableLayerReadsARun)
{
    alps::hdf5::archive archive(fixtures + "legacy_observables.h5", "r");
    const alps::alea::result_set results = alps::alea::load_results(archive, "/run");
    EXPECT_FALSE(results.has("Sign * O"));

    numbers got;
    const alps::alea::observable_result x = results["x"], s = results["s"], v = results["v"], o = results["O"];
    got["x/count"] = x.count();
    got["x/mean"] = x.mean();
    got["x/error"] = x.error();
    got["x/variance"] = x.variance_vector()[0];
    got["x/tau"] = x.tau();
    got["x/converged"] = x.converged();
    got["x/bins"] = x.bin_means().size();
    got["x/binsize"] = x.bin_size();
    // The file of a running 3.0 simulation holds bin sums.
    got["x/bin_first"] = x.bin_means().front()[0] * 256;
    got["x/bin_last"] = x.bin_means().back()[0] * 256;
    got["s/count"] = s.count();
    got["s/mean"] = s.mean();
    got["s/error"] = s.error();
    got["v/count"] = v.count();
    for (int i = 0; i < 2; ++i) {
        got["v/mean/" + std::to_string(i)] = v.mean_vector()[i];
        got["v/error/" + std::to_string(i)] = v.error_vector()[i];
        got["v/tau/" + std::to_string(i)] = v.tau_vector()[i];
    }
    got["Sign/mean"] = results["Sign"].mean();
    // ALEA reports the samples in full bins, 78 * 256, as the count of a
    // signed observable. The file has the number of samples.
    EXPECT_EQ(o.count(), 20000u);
    got["O/mean"] = o.mean();
    got["O/error"] = o.error();
    got["h/count"] = results["h"].count();
    for (int bin = 0; bin < 10; ++bin)
        got["h/" + std::to_string(bin)] = results["h"].histogram()[bin];

    numbers want = expected("run");
    want.erase("x/error_level0");
    want.erase("x/error_level3");
    want.erase("O/count");
    expect_equal(want, got);

    // The stored bins serve functions of results, also for the signed one.
    EXPECT_FALSE(s.has_bins());
    EXPECT_NEAR((x * x).mean(), x.mean() * x.mean(), 1e-3);
    EXPECT_NEAR((o * 2.0).error(), 2 * o.error(), 0.02 * o.error());

    // Measurements of an earlier version cannot be continued, and say so.
    alps::alea::observable_set set;
    EXPECT_THROW(alps::alea::load_checkpoint(archive, "/run", set), std::runtime_error);
}

// Merged runs, as a task's result group holds them, read back with the stored
// numbers. The group with a signed observable, which the ALEA of 3.0 wrote but
// could not read, reads too.
TEST(FrozenFiles, ObservableLayerReadsMergedRuns)
{
    alps::hdf5::archive archive(fixtures + "legacy_observables.h5", "r");
    const alps::alea::result_set results = alps::alea::load_results(archive, "/merged");
    numbers got;
    const alps::alea::observable_result x = results["x"], v = results["v"], h = results["h"];
    got["x/count"] = x.count();
    got["x/mean"] = x.mean();
    got["x/error"] = x.error();
    got["x/tau"] = x.tau();
    got["x/bins"] = x.bin_means().size();
    for (int i = 0; i < 2; ++i) {
        got["v/mean/" + std::to_string(i)] = v.mean_vector()[i];
        got["v/error/" + std::to_string(i)] = v.error_vector()[i];
    }
    got["h/count"] = h.count();
    got["h/0"] = h.histogram()[0];
    got["h/9"] = h.histogram()[9];
    expect_equal(expected("merged"), got);

    const alps::alea::result_set signed_runs = alps::alea::load_results(archive, "/merged_signed");
    EXPECT_EQ(signed_runs["O"].sign(), "Sign");
    // Signed average 1.5; two runs of 20000 samples give an error near 0.009.
    EXPECT_NEAR(signed_runs["O"].mean(), 1.5, 0.05);
}
