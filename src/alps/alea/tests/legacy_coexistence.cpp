// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// Legacy ALEA and ALPSCore's ALEA share the include prefix alps/alea/ and the
// namespace alps::alea. Check that one program can use both, and that their
// results can live side by side in one archive.
#include <alps/alea.h>
#include <alps/alea/mcanalyze.hpp>
#include <alps/alea/mcdata.hpp>
#include <alps/alea/value_with_error.hpp>

#include <alps/alea.hpp>

#include <alps/hdf5.hpp>
#include <alps/testing/temporary_directory.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>
#include <random>
#include <vector>

namespace {

std::vector<double> samples()
{
    std::mt19937 engine(42);
    std::normal_distribution<double> normal(1.0, 0.5);
    std::vector<double> values(4096);
    for (double & value : values)
        value = normal(engine);
    return values;
}

}

TEST(AleaLegacyCoexistence, OneArchiveHoldsBoth)
{
    const std::vector<double> values = samples();
    const double count = static_cast<double>(values.size());
    const double mean = std::accumulate(values.begin(), values.end(), 0.) / count;

    alps::ObservableSet legacy;
    legacy << alps::RealObservable("x");
    alps::alea::var_acc<double> variance(1);
    alps::alea::batch_acc<double> batch(1);
    for (double value : values) {
        legacy["x"] << value;
        variance << value;
        batch << value;
    }
    const alps::alea::var_result<double> variance_result = variance.finalize();
    const alps::alea::batch_result<double> batch_result = batch.finalize();

    alps::testing::TemporaryDirectory temporary;
    const std::string filename = (temporary.path() / "coexistence.h5").string();
    {
        alps::hdf5::archive archive(filename, "w");
        archive["/legacy"] << legacy;
        alps::alea::hdf5_serializer serializer(archive, "/modern");
        alps::alea::serialize(serializer, "variance", variance_result);
        alps::alea::serialize(serializer, "batch", batch_result);
    }

    alps::ObservableSet legacy_loaded;
    alps::alea::var_result<double> variance_loaded;
    alps::alea::batch_result<double> batch_loaded;
    {
        alps::hdf5::archive archive(filename, "r");
        EXPECT_TRUE(archive.is_group("/legacy"));
        EXPECT_TRUE(archive.is_group("/modern/variance"));
        EXPECT_TRUE(archive.is_group("/modern/batch"));
        archive["/legacy"] >> legacy_loaded;
        alps::alea::hdf5_serializer serializer(archive, "/modern");
        alps::alea::deserialize(serializer, "variance", variance_loaded);
        alps::alea::deserialize(serializer, "batch", batch_loaded);
    }

    const alps::RealObservable & x = legacy_loaded.get<alps::RealObservable>("x");
    EXPECT_EQ(x.count(), values.size());
    EXPECT_EQ(variance_loaded.count(), values.size());
    EXPECT_EQ(batch_loaded.count(), values.size());
    EXPECT_NEAR(x.mean(), mean, 1e-12);
    EXPECT_NEAR(variance_loaded.mean()[0], mean, 1e-12);
    EXPECT_NEAR(batch_loaded.mean()[0], mean, 1e-12);
    // Both report the naive standard error; they may differ by the n / (n - 1)
    // convention, a relative 1e-4 for this sample.
    EXPECT_NEAR(variance_loaded.stderror()[0], x.error(0), 1e-3 * x.error(0));
}

TEST(AleaLegacyCoexistence, LegacyAnalysisNamesStillResolve)
{
    // Legacy alps::alea functions and value types, next to the modern classes.
    const std::vector<double> values = samples();
    alps::RealObservable legacy("x");
    alps::alea::mean_acc<double> modern(1);
    for (double value : values) {
        legacy << value;
        modern << value;
    }
    const alps::alea::mcdata<double> data(legacy);
    const alps::alea::value_with_error<double> estimate(data.mean(), data.error());
    EXPECT_EQ(estimate.mean(), data.mean());
    EXPECT_GT(estimate.error(), 0.);
    EXPECT_NEAR(modern.finalize().mean()[0], data.mean(), 1e-12);
}
