// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
// A sign-problem measurement with a known answer. The sign s is +1 with
// probability 0.8 and -1 otherwise (mean 0.6); O = 1 + 0.3 s + z with z a
// standard normal deviate. The simulation accumulates s O and s, and the
// signed average is <s O> / <s> = 0.9 / 0.6 = 1.5. The delta method with
// Var(s O) = 1.64, Var(s) = 0.64, Cov(s O, s) = 0.64 gives a variance of
// 3.2222 / N. Treating numerator and sign as independent would give an error
// 1.63 times too large.
#include "streams.hpp"

#include <alps/alea.h>
#include <alps/hdf5.hpp>
#include <alps/testing/temporary_directory.hpp>

#include <gtest/gtest.h>

#include <cmath>

namespace {

const std::size_t n = std::size_t(1) << 16;
const double exact = 1.5;
const double exact_error = std::sqrt((1.64 / 0.36 + 0.81 * 0.64 / 0.1296 - 2 * 0.9 * 0.64 / 0.216) / n);

alps::ObservableSet measure() {
    alps::ObservableSet set;
    set << alps::RealObservable("Sign") << alps::SignedObservable<alps::RealObservable>("O");
    statistics_test::normal xi(8128);
    for (std::size_t t = 0; t < n; ++t) {
        const double sign = xi.uniform() < 0.8 ? 1. : -1.;
        set["Sign"] << sign;
        set["O"] << sign * (1 + 0.3 * sign + xi());
    }
    return set;
}

}

TEST(SignedRatio, Alea)
{
    alps::ObservableSet set = measure();
    const alps::RealObsevaluator ratio = set["O"];
    EXPECT_NEAR(exact, ratio.mean(), 4 * exact_error);
    // A jackknife over 128 bins scatters by 6%.
    EXPECT_NEAR(1.0, ratio.error() / exact_error, 0.25);

    // The value, its error and the name of its sign survive a file.
    alps::testing::TemporaryDirectory temporary;
    const std::string filename = (temporary.path() / "signed.h5").string();
    {
        alps::hdf5::archive archive(filename, "w");
        archive["/set"] << set;
    }
    alps::ObservableSet loaded;
    {
        alps::hdf5::archive archive(filename, "r");
        archive["/set"] >> loaded;
    }
    const alps::RealObsevaluator reloaded = loaded["O"];
    EXPECT_EQ(ratio.mean(), reloaded.mean());
    EXPECT_EQ(ratio.error(), reloaded.error());
    EXPECT_EQ(loaded.get<alps::SignedObservable<alps::RealObservable> >("O").sign_name(), "Sign");
}
