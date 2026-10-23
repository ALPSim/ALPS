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

TEST(FulldiagParameters, FieldCanBeAnExpression) {
    // T, beta, H and h are read with the model parameters when the XML output
    // is written, so h = J/2 does not stop the task
    alps::params const p = read("T = 0.5\nh = J/2");
    EXPECT_FALSE(p.exists("T") || p.exists("h"));
    alps::params const b = read("beta = 2*J\nH = J/2");
    EXPECT_FALSE(b.exists("beta") || b.exists("H"));
}

TEST(FulldiagParameters, NoLanczosIterations) {
    EXPECT_FALSE(read("MAX_ITERATIONS = abc").exists("MAX_ITERATIONS"));
}

}
