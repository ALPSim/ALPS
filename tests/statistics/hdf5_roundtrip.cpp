// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Saving and loading. HDF5 stores doubles exactly, so a reloaded object must
// report exactly what the saved one did, and a run that is saved, reloaded and
// continued must end exactly where an uninterrupted run ends.
#include "streams.hpp"
#include "fixtures/legacy_observables.hpp"

#include <alps/alea.h>
#include <alps/alea.hpp>
#include <alps/accumulators.hpp>
#include <alps/hdf5.hpp>
#include <alps/testing/temporary_directory.hpp>

#include <gtest/gtest.h>

#include <string>
#include <vector>

namespace {

using statistics_test::describe_run;
using statistics_test::numbers;

// Deliberately not a power of two: every level holds an unfinished bin.
const std::size_t first_part = 5000, total = 8000;

void expect_equal(const numbers & want, const numbers & got) {
    ASSERT_EQ(want.size(), got.size());
    for (const auto & value : want) {
        if (value.second != value.second)
            EXPECT_NE(got.at(value.first), got.at(value.first)) << value.first;
        else
            EXPECT_EQ(value.second, got.at(value.first)) << value.first;
    }
}

template <typename Set> Set reload(const Set & set, const std::string & filename) {
    {
        alps::hdf5::archive archive(filename, "w");
        archive["/set"] << set;
    }
    Set loaded;
    alps::hdf5::archive archive(filename, "r");
    archive["/set"] >> loaded;
    return loaded;
}

}

TEST(Hdf5RoundTrip, Alea)
{
    alps::testing::TemporaryDirectory temporary;
    const alps::ObservableSet set = statistics_test::legacy_run(9, total);
    expect_equal(describe_run(set), describe_run(reload(set, (temporary.path() / "alea.h5").string())));
}

TEST(Hdf5RoundTrip, AleaResume)
{
    alps::testing::TemporaryDirectory temporary;
    const std::vector<double> values = statistics_test::ar1_series(9, 0.7, total);
    alps::ObservableSet whole, part;
    whole << alps::RealObservable("x") << alps::SimpleRealObservable("s");
    part << alps::RealObservable("x") << alps::SimpleRealObservable("s");
    for (std::size_t t = 0; t < total; ++t) {
        whole["x"] << values[t];
        whole["s"] << values[t];
        if (t < first_part) {
            part["x"] << values[t];
            part["s"] << values[t];
        }
    }
    alps::ObservableSet resumed = reload(part, (temporary.path() / "alea.h5").string());
    for (std::size_t t = first_part; t < total; ++t) {
        resumed["x"] << values[t];
        resumed["s"] << values[t];
    }
    const alps::RealObservable & a = whole.get<alps::RealObservable>("x");
    const alps::RealObservable & b = resumed.get<alps::RealObservable>("x");
    EXPECT_EQ(a.count(), b.count());
    EXPECT_EQ(a.mean(), b.mean());
    EXPECT_EQ(a.error(), b.error());
    EXPECT_EQ(a.tau(), b.tau());
    for (unsigned level = 0; level < 6; ++level)
        EXPECT_EQ(a.error(level), b.error(level)) << "level " << level;
    ASSERT_EQ(a.bin_number(), b.bin_number());
    for (unsigned bin = 0; bin < a.bin_number(); ++bin)
        EXPECT_EQ(a.bin_value(bin), b.bin_value(bin)) << "bin " << bin;
    EXPECT_EQ(whole.get<alps::SimpleRealObservable>("s").error(), resumed.get<alps::SimpleRealObservable>("s").error());
}

TEST(Hdf5RoundTrip, AccumulatorsResume)
{
    namespace acc = alps::accumulators;
    alps::testing::TemporaryDirectory temporary;
    const std::vector<double> values = statistics_test::ar1_series(9, 0.7, total);
    acc::accumulator_set whole, part;
    for (acc::accumulator_set * set : {&whole, &part})
        *set << acc::NoBinningAccumulator<double>("simple") << acc::LogBinningAccumulator<double>("log")
             << acc::FullBinningAccumulator<double>("full");
    const char * names[] = {"simple", "log", "full"};
    for (std::size_t t = 0; t < total; ++t)
        for (const char * name : names) {
            whole[name] << values[t];
            if (t < first_part)
                part[name] << values[t];
        }
    acc::accumulator_set resumed = reload(part, (temporary.path() / "accumulators.h5").string());
    for (std::size_t t = first_part; t < total; ++t)
        for (const char * name : names)
            resumed[name] << values[t];
    const acc::result_set a(whole), b(resumed);
    for (const char * name : names) {
        EXPECT_EQ(a[name].count(), b[name].count()) << name;
        EXPECT_EQ(a[name].mean<double>(), b[name].mean<double>()) << name;
        EXPECT_EQ(a[name].error<double>(), b[name].error<double>()) << name;
    }
    EXPECT_EQ(a["log"].autocorrelation<double>(), b["log"].autocorrelation<double>());
}

TEST(Hdf5RoundTrip, AleaV2Results)
{
    using namespace alps::alea;
    alps::testing::TemporaryDirectory temporary;
    const std::string filename = (temporary.path() / "aleav2.h5").string();
    var_acc<double> var(1);
    autocorr_acc<double> autocorr(1);
    batch_acc<double> batch(1);
    cov_acc<double> cov(2);
    statistics_test::ar1 source(9, 0.7);
    for (std::size_t t = 0; t < total; ++t) {
        const double x = source();
        var << x;
        autocorr << x;
        batch << x;
        cov << std::vector<double>({x, x * x});
    }
    const var_result<double> var_saved = var.finalize();
    const autocorr_result<double> autocorr_saved = autocorr.finalize();
    const batch_result<double> batch_saved = batch.finalize();
    const cov_result<double> cov_saved = cov.finalize();
    {
        alps::hdf5::archive archive(filename, "w");
        hdf5_serializer serializer(archive, "/results");
        serialize(serializer, "var", var_saved);
        serialize(serializer, "autocorr", autocorr_saved);
        serialize(serializer, "batch", batch_saved);
        serialize(serializer, "cov", cov_saved);
    }
    var_result<double> var_loaded;
    autocorr_result<double> autocorr_loaded;
    batch_result<double> batch_loaded;
    cov_result<double> cov_loaded;
    {
        alps::hdf5::archive archive(filename, "r");
        hdf5_serializer serializer(archive, "/results");
        deserialize(serializer, "var", var_loaded);
        deserialize(serializer, "autocorr", autocorr_loaded);
        deserialize(serializer, "batch", batch_loaded);
        deserialize(serializer, "cov", cov_loaded);
    }
    EXPECT_EQ(var_saved.count(), var_loaded.count());
    EXPECT_EQ(var_saved.mean()[0], var_loaded.mean()[0]);
    EXPECT_EQ(var_saved.stderror()[0], var_loaded.stderror()[0]);
    EXPECT_EQ(autocorr_saved.count(), autocorr_loaded.count());
    EXPECT_EQ(autocorr_saved.mean()[0], autocorr_loaded.mean()[0]);
    EXPECT_EQ(autocorr_saved.stderror()[0], autocorr_loaded.stderror()[0]);
    EXPECT_EQ(autocorr_saved.tau()[0], autocorr_loaded.tau()[0]);
    ASSERT_EQ(autocorr_saved.nlevel(), autocorr_loaded.nlevel());
    for (std::size_t level = 0; level < autocorr_saved.nlevel(); ++level)
        EXPECT_EQ(autocorr_saved.level(level).stderror()[0], autocorr_loaded.level(level).stderror()[0]) << level;
    EXPECT_EQ(batch_saved.count(), batch_loaded.count());
    EXPECT_EQ(batch_saved.mean()[0], batch_loaded.mean()[0]);
    EXPECT_EQ(batch_saved.stderror()[0], batch_loaded.stderror()[0]);
    EXPECT_EQ(cov_saved.cov()(0, 1), cov_loaded.cov()(0, 1));
}

// The five ALEA v2 accumulators, saved in the middle of a run and continued
// from the file, end bit for bit where the uninterrupted ones do.
TEST(Hdf5RoundTrip, AleaV2Resume)
{
    using namespace alps::alea;
    alps::testing::TemporaryDirectory temporary;
    const std::string filename = (temporary.path() / "aleav2_state.h5").string();
    mean_acc<double> mean(1), mean_resumed;
    var_acc<double> var(1, 3), var_resumed;
    cov_acc<double> cov(2), cov_resumed;
    autocorr_acc<double> autocorr(1), autocorr_resumed;
    batch_acc<double> batch(1, 16), batch_resumed;
    statistics_test::ar1 source(9, 0.7);
    for (std::size_t t = 0; t < total; ++t) {
        if (t == first_part) {
            {
                alps::hdf5::archive archive(filename, "w");
                hdf5_serializer serializer(archive, "/checkpoint");
                mean.save_state(serializer, "mean");
                var.save_state(serializer, "var");
                cov.save_state(serializer, "cov");
                autocorr.save_state(serializer, "autocorr");
                batch.save_state(serializer, "batch");
            }
            alps::hdf5::archive archive(filename, "r");
            hdf5_serializer serializer(archive, "/checkpoint");
            mean_resumed.load_state(serializer, "mean");
            var_resumed.load_state(serializer, "var");
            cov_resumed.load_state(serializer, "cov");
            autocorr_resumed.load_state(serializer, "autocorr");
            batch_resumed.load_state(serializer, "batch");
            // A state is not mistaken for one of another accumulator.
            EXPECT_THROW(var_resumed.load_state(serializer, "batch"), std::runtime_error);
        }
        const double x = source();
        const std::vector<double> pair = {x, x * x};
        mean << x;
        var << x;
        cov << pair;
        autocorr << x;
        batch << x;
        if (t >= first_part) {
            mean_resumed << x;
            var_resumed << x;
            cov_resumed << pair;
            autocorr_resumed << x;
            batch_resumed << x;
        }
    }
    EXPECT_TRUE(mean.finalize() == mean_resumed.finalize());
    EXPECT_TRUE(var.finalize() == var_resumed.finalize());
    EXPECT_TRUE(cov.finalize() == cov_resumed.finalize());
    EXPECT_TRUE(autocorr.finalize() == autocorr_resumed.finalize());
    const batch_result<double> batch_result_ = batch.finalize(), batch_resumed_result = batch_resumed.finalize();
    EXPECT_TRUE(batch_result_ == batch_resumed_result);
    EXPECT_EQ(batch_result_.store().batch_size(), batch_resumed_result.store().batch_size());
    EXPECT_EQ(batch_result_.count(), total);
}

// Samples added directly, without a batch, give the same sums as samples added
// through an adapter, and their state (an empty batch is not written) resumes.
TEST(Hdf5RoundTrip, AleaDirectSamples)
{
    using namespace alps::alea;
    alps::testing::TemporaryDirectory temporary;
    const std::string filename = (temporary.path() / "alea_direct.h5").string();
    var_acc<double> direct(2), adapted(2), resumed, batched_direct(2, 3), batched(2, 3);
    statistics_test::ar1 source(11, 0.5);
    for (std::size_t t = 0; t < total; ++t) {
        if (t == first_part) {
            {
                alps::hdf5::archive archive(filename, "w");
                hdf5_serializer serializer(archive, "/checkpoint");
                direct.save_state(serializer, "var");
            }
            alps::hdf5::archive archive(filename, "r");
            hdf5_serializer serializer(archive, "/checkpoint");
            resumed.load_state(serializer, "var");
        }
        const double x = source();
        const std::vector<double> pair = {x, -2 * x};
        direct.add_sample(pair.data());
        adapted << make_adapter(pair);
        batched_direct.add_sample(pair.data());
        batched << make_adapter(pair);
        if (t >= first_part)
            resumed.add_sample(pair.data());
    }
    EXPECT_TRUE(direct.result() == adapted.result());
    EXPECT_TRUE(direct.result() == resumed.result());
    EXPECT_TRUE(batched_direct.result() == batched.result());
    EXPECT_EQ(direct.count(), total);
}
