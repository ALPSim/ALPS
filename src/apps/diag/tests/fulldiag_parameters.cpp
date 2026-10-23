// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include "../fulldiag/fulldiag.h"

#include <alps/parameter.h>
#include <gtest/gtest.h>

#include <sstream>
#include <string>

namespace {

typedef FullDiagMatrix<double> matrix_type;

alps::params read(std::string const & text) {
    std::istringstream in(text);
    alps::Parameters parameters;
    in >> parameters;
    return matrix_type::read_parameters(parameters, &matrix_type::define_parameters);
}

TEST(FulldiagParameters, FieldOnlyWithTemperature) {
    alps::params const none = read("h = J/2");
    EXPECT_FALSE(none.exists("T") || none.exists("beta") || none.exists("h"));

    alps::params const t = read("T = 0.5\nbeta = 3\nh = 2");
    EXPECT_EQ(t["T"].as<double>(), 0.5);
    EXPECT_FALSE(t.exists("beta"));
    EXPECT_EQ(t["h"].as<double>(), 2.);

    alps::params const b = read("beta = 3\nH = 1\nh = J/2");
    EXPECT_EQ(b["beta"].as<double>(), 3.);
    EXPECT_EQ(b["H"].as<double>(), 1.);
    EXPECT_FALSE(b.exists("h"));
}

TEST(FulldiagParameters, NoLanczosIterations) {
    EXPECT_FALSE(read("MAX_ITERATIONS = abc").exists("MAX_ITERATIONS"));
}

}
