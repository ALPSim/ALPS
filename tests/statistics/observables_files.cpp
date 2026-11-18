// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// The observable layer in files and messages: results written and read back,
// checkpoints, and results sent to another process.
#include "streams.hpp"

#include <alps/alea/observables.hpp>
#include <alps/alea/observables_io.hpp>
#include <alps/hdf5.hpp>
#include <alps/testing/temporary_directory.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace {

using namespace alps::alea;

observable_set declare() {
    observable_set set;
    set << observable("x") << observable("plain").without_bins() << observable("v").labels({"a", "b"})
        << observable("Sign") << observable("O").signed_by("Sign") << histogram("h", 0, 9) << observable("unused");
    return set;
}

// Samples first to last - 1 of one fixed stream.
void measure(observable_set & set, std::size_t first, std::size_t last) {
    statistics_test::ar1 correlated(77, 0.7);
    statistics_test::normal xi(78);
    for (std::size_t t = 0; t < last; ++t) {
        const double x = correlated();
        const double sign = xi.uniform() < 0.8 ? 1. : -1.;
        const double o = sign * (1 + 0.3 * sign + xi());
        const int bin = static_cast<int>(10 * xi.uniform());
        if (t < first)
            continue;
        set["x"] << x;
        set["plain"] << x;
        set["v"] << std::vector<double>({x, 2 * x + o});
        set["Sign"] << sign;
        set["O"] << o;
        set["h"] << bin;
    }
}

void expect_same(const std::vector<double> & a, const std::vector<double> & b) {
    ASSERT_EQ(a.size(), b.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] != a[i])
            EXPECT_NE(b[i], b[i]);
        else
            EXPECT_EQ(a[i], b[i]);
    }
}

// Two sets of results that must agree in every number, bit for bit.
void expect_same(const result_set & a, const result_set & b) {
    ASSERT_EQ(a.names(), b.names());
    for (const std::string & name : a.names()) {
        SCOPED_TRACE(name);
        ASSERT_EQ(a[name].count(), b[name].count());
        ASSERT_EQ(a[name].size(), b[name].size());
        if (name == "h") {
            EXPECT_EQ(a[name].histogram(), b[name].histogram());
            continue;
        }
        if (a[name].count() == 0)
            continue;
        expect_same(a[name].mean_vector(), b[name].mean_vector());
        expect_same(a[name].error_vector(), b[name].error_vector());
        expect_same(a[name].tau_vector(), b[name].tau_vector());
        expect_same(a[name].variance_vector(), b[name].variance_vector());
        EXPECT_EQ(a[name].converged_vector(), b[name].converged_vector());
        EXPECT_EQ(a[name].labels(), b[name].labels());
        ASSERT_EQ(a[name].has_bins(), b[name].has_bins());
        expect_same(a[name].binning_errors(), b[name].binning_errors());
    }
}

}

// 2^13 samples fill 128 bins of 64 exactly, so the file holds every sample
// and a function of the results read back is the function of the originals.
TEST(ObservableFiles, ResultsReadBackAsWritten)
{
    alps::testing::TemporaryDirectory temporary;
    const std::string filename = (temporary.path() / "results.h5").string();
    observable_set set = declare();
    measure(set, 0, std::size_t(1) << 13);
    const result_set written = set.results();
    {
        alps::hdf5::archive archive(filename, "w");
        save_results(archive, "/simulation/results", written);
    }
    alps::hdf5::archive archive(filename, "r");
    const result_set read = load_results(archive, "/simulation/results");
    expect_same(written, read);
    EXPECT_EQ(read["O"].sign(), "Sign");
    EXPECT_EQ(read["x"].bin_size(), 64u);
    ASSERT_EQ(read["x"].bin_means().size(), 128u);
    EXPECT_NEAR(read["x"].bin_means()[5][0], written["x"].bin_means()[5][0], 1e-15);

    const observable_result ratio = written["x"] / written["O"], ratio_read = read["x"] / read["O"];
    EXPECT_NEAR(ratio.mean(), ratio_read.mean(), 1e-12);
    EXPECT_NEAR(ratio.error(), ratio_read.error(), 1e-12);
    EXPECT_THROW(read["x"] * read["plain"], std::invalid_argument);

    // The layout of SCHEMA.md: what a 3.0 reader looks for is where it was,
    // and the internals of the old accumulators are gone.
    for (const char * path : {"count", "mean/value", "mean/error", "mean/error_convergence", "variance/value",
                              "tau/value", "timeseries/data", "binning/binsize", "binning/error"})
        EXPECT_TRUE(archive.is_data(std::string("/simulation/results/x/") + path)) << path;
    for (const char * path : {"sum", "sum2", "timeseries/logbinning", "timeseries/data2", "timeseries/partialbin"})
        EXPECT_FALSE(archive.is_data(std::string("/simulation/results/x/") + path)) << path;
    EXPECT_TRUE(archive.is_scalar("/simulation/results/x/mean/value"));
    EXPECT_FALSE(archive.is_scalar("/simulation/results/v/mean/value"));
    EXPECT_TRUE(archive.is_attribute("/simulation/results/x/timeseries/data/@binsize"));
    EXPECT_TRUE(archive.is_data("/simulation/results/Sign * O/timeseries/data"));
    EXPECT_TRUE(archive.is_attribute("/simulation/results/O/@sign"));
    EXPECT_TRUE(archive.is_data("/simulation/results/h/histogram"));

    // Stored results are not accumulated data: no merging, no continuing.
    EXPECT_THROW(set.merge(read), std::invalid_argument);
    observable_set other;
    EXPECT_THROW(load_checkpoint(archive, "/simulation/results", other), std::runtime_error);
}

// Saved after 5000 of 8000 samples and continued from the file, a run ends
// bit for bit where the uninterrupted run ends.
TEST(ObservableFiles, CheckpointContinuesExactly)
{
    alps::testing::TemporaryDirectory temporary;
    const std::string filename = (temporary.path() / "checkpoint.h5").string();
    observable_set whole = declare(), first = declare();
    measure(whole, 0, 8000);
    measure(first, 0, 5000);
    {
        alps::hdf5::archive archive(filename, "w");
        save_checkpoint(archive, "/checkpoint/measurements", first);
    }
    observable_set resumed;
    resumed << observable("x") << observable("extra");
    {
        alps::hdf5::archive archive(filename, "r");
        load_checkpoint(archive, "/checkpoint/measurements", resumed);
    }
    EXPECT_TRUE(resumed.has("extra"));
    EXPECT_EQ(resumed["x"].count(), 5000u);
    measure(resumed, 5000, 8000);
    whole << observable("extra");
    expect_same(whole.results(), resumed.results());
    const observable_result a = whole.results()["x"] * whole.results()["O"],
                            b = resumed.results()["x"] * resumed.results()["O"];
    EXPECT_EQ(a.mean(), b.mean());
    EXPECT_EQ(a.error(), b.error());
}

// Results sent as bytes arrive as they were and can be merged on arrival.
TEST(ObservableFiles, PackedResultsMergeLikeTheOriginals)
{
    observable_set first = declare(), second = declare();
    measure(first, 0, 3000);
    measure(second, 3000, 8000);
    const result_set sent = second.results();
    const result_set received = unpack(pack(sent));
    expect_same(sent, received);

    observable_set direct = first, by_message = first;
    direct.merge(sent);
    by_message.merge(received);
    expect_same(direct.results(), by_message.results());
    EXPECT_EQ(direct.results()["x"].count(), 8000u);
    EXPECT_THROW(unpack("not packed results"), std::runtime_error);
}

// Bins of earlier versions: an evaluated result (marked "@cannotrebin", as
// ALEA's mcdata and ALPSCore write it) holds bin means, the observable of a
// running simulation bin sums.
TEST(ObservableFiles, EarlierBinsAreMeansOrSums)
{
    alps::testing::TemporaryDirectory temporary;
    const std::string name = (temporary.path() / "earlier.h5").string();
    {
        alps::hdf5::archive archive(name, "w");
        for (const char * group : {"/r/evaluated", "/r/running"}) {
            const std::string path = group;
            const bool evaluated = path == "/r/evaluated";
            archive[path + "/count"] << std::uint64_t(64);
            archive[path + "/mean/value"] << 2.5;
            archive[path + "/mean/error"] << 0.5;
            archive[path + "/timeseries/data"] << (evaluated ? std::vector<double>({1., 2., 3., 4.})
                                                             : std::vector<double>({16., 32., 48., 64.}));
            archive[path + "/timeseries/data/@binsize"] << std::uint64_t(16);
            if (evaluated)
                archive[path + "/@cannotrebin"] << false;
        }
    }
    alps::hdf5::archive archive(name, "r");
    const result_set results = load_results(archive, "/r");
    for (const char * group : {"evaluated", "running"}) {
        SCOPED_TRACE(group);
        const observable_result r = results[group];
        ASSERT_TRUE(r.has_bins());
        EXPECT_EQ(r.bin_size(), 16u);
        const std::vector<std::vector<double> > bins = r.bin_means();
        ASSERT_EQ(bins.size(), 4u);
        for (std::size_t i = 0; i < 4; ++i)
            EXPECT_EQ(bins[i][0], i + 1.);
    }
}
