// Copyright (C) 2026 by the ALPS collaboration
// SPDX-License-Identifier: MIT
#include <alps/input/ctint.hpp>
#include <alps/input/hybridization.hpp>
#include <alps/parameter_input.hpp>
#include <gtest/gtest.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

// CT-INT parameters as the solver gets them: supplied text, then defined.
alps::input::ctint_parameters ctint(std::string const & supplied) {
    std::istringstream in(supplied);
    alps::params parameters;
    for (alps::Parameter const & p : alps::read_parameter_text(in))
        parameters.supply(p.key(), p.value().c_str());
    alps::input::define_ctint_parameters(parameters);
    std::ostringstream missing;
    if (parameters.has_missing(missing))
        throw std::invalid_argument(missing.str());
    return alps::input::read_ctint_parameters(parameters);
}

std::string const ctint_required =
    "SWEEPS = 1000; THERMALIZATION = 10; MAX_TIME = 60; BETA = 4; U = 2; ALPHA = -0.01; INFILE = \"g0.h5\"\n";

TEST(CtintInput, Defaults) {
    alps::input::ctint_parameters const p = ctint(ctint_required + "N_MATSUBARA = 40; N_TAU = 50\n");
    EXPECT_EQ(p.max_order, 2048u);
    EXPECT_EQ(p.flavors, 2u);
    EXPECT_EQ(p.sites, 1u);
    EXPECT_EQ(p.n_matsubara, 40u);
    EXPECT_EQ(p.n_matsubara_measurements, 40u);
    EXPECT_EQ(p.n_tau, 50u);
    EXPECT_EQ(p.n_self, 500u);
    EXPECT_EQ(p.sweeps, 1000u);
    EXPECT_EQ(p.thermalization, 10u);
    EXPECT_EQ(p.max_time, 60);
    EXPECT_EQ(p.recalc_period, 5000u);
    EXPECT_EQ(p.measurement_period, 200u);
    EXPECT_EQ(p.convergence_check_period, 5000u);
    EXPECT_EQ(p.almost_zero, 1e-16);
    EXPECT_FALSE(p.histogram_measurement);
    EXPECT_FALSE(p.atomic);
    EXPECT_FALSE(p.dmft_framework);
    EXPECT_EQ(p.tau_discretization_for_exp, 0u);
    EXPECT_EQ(p.infile, "g0.h5");
}

// NMATSUBARA and N, the older names, take precedence; FLAVORS over N_ORBITALS.
TEST(CtintInput, OlderNamesTakePrecedence) {
    alps::input::ctint_parameters const p = ctint(ctint_required
        + "NMATSUBARA = 30; N_MATSUBARA = 40; N = 20; N_TAU = 50; FLAVORS = 4; N_ORBITALS = 3\n");
    EXPECT_EQ(p.n_matsubara, 30u);
    EXPECT_EQ(p.n_tau, 20u);
    EXPECT_EQ(p.n_self, 200u);
    EXPECT_EQ(p.flavors, 4u);
}

TEST(CtintInput, PeriodsAndSizesFollowTheirAlternatives) {
    alps::input::ctint_parameters const p = ctint(ctint_required
        + "N_MATSUBARA = 40; N_TAU = 50; N_ORBITALS = 3; N_MEAS = 7; RECALC_PERIOD = 333\n");
    EXPECT_EQ(p.flavors, 3u);
    EXPECT_EQ(p.measurement_period, 7u);
    EXPECT_EQ(p.convergence_check_period, 333u);
    alps::input::ctint_parameters const q = ctint(ctint_required
        + "N_MATSUBARA = 40; N_TAU = 50; NMATSUBARA_MEASUREMENTS = 12; NSELF = 123;"
          " MEASUREMENT_PERIOD = 17; CONVERGENCE_CHECK_PERIOD = 999; TAU_DISCRETIZATION_FOR_EXP = 64\n");
    EXPECT_EQ(q.n_matsubara_measurements, 12u);
    EXPECT_EQ(q.n_self, 123u);
    EXPECT_EQ(q.measurement_period, 17u);
    EXPECT_EQ(q.convergence_check_period, 999u);
    EXPECT_EQ(q.tau_discretization_for_exp, 64u);
}

TEST(CtintInput, InfileIsRequiredUnlessAtomic) {
    std::string const sizes = "SWEEPS = 1000; THERMALIZATION = 10; MAX_TIME = 60; BETA = 4; U = 2; ALPHA = 0;"
                              " N_MATSUBARA = 40; N_TAU = 50\n";
    EXPECT_THROW(ctint(sizes), std::invalid_argument);
    alps::input::ctint_parameters const p = ctint(sizes + "ATOMIC = 1\n");
    EXPECT_TRUE(p.atomic);
    EXPECT_EQ(p.infile, "");
}

// CT-HYB parameters as the solver gets them: supplied text, then defined.
alps::input::hybridization_parameters cthyb(std::string const & supplied) {
    std::istringstream in(supplied);
    alps::params parameters;
    for (alps::Parameter const & p : alps::read_parameter_text(in))
        parameters.supply(p.key(), p.value().c_str());
    alps::input::define_hybridization_parameters(parameters);
    std::ostringstream missing;
    if (parameters.has_missing(missing))
        throw std::invalid_argument(missing.str());
    return alps::input::read_hybridization_parameters(parameters);
}

std::string const cthyb_required =
    "SWEEPS = 1000; THERMALIZATION = 10; MAX_TIME = 60; N_MEAS = 100; BETA = 4;"
    " N_TAU = 50; N_ORBITALS = 2; DELTA = \"delta.dat\"\n";

TEST(HybridizationInput, Defaults) {
    alps::input::hybridization_parameters const p = cthyb(cthyb_required + "U = 3\n");
    EXPECT_EQ(p.output_period, 100000);
    EXPECT_EQ(p.n_histogram_orders, 50);
    EXPECT_TRUE(p.measure_time);
    EXPECT_FALSE(p.measure_freq);
    EXPECT_TRUE(p.accurate_covariance);
    EXPECT_EQ(p.basepath, "");
    EXPECT_EQ(p.mu, 0.);
    EXPECT_EQ(p.J, 0.);
    EXPECT_EQ(p.U, 3.);
    EXPECT_EQ(p.U_prime, 3.);
    EXPECT_FALSE(p.u_ignored);
    EXPECT_FALSE(p.mu_ignored);
}

TEST(HybridizationInput, UPrimeDefaultsToUMinusTwoJ) {
    EXPECT_EQ(cthyb(cthyb_required + "U = 3; J = 0.5\n").U_prime, 2.);
    EXPECT_EQ(cthyb(cthyb_required + "U = 3; J = 0.5; U' = 1.75\n").U_prime, 1.75);
}

TEST(HybridizationInput, UMatrixReplacesU) {
    EXPECT_THROW(cthyb(cthyb_required), std::invalid_argument);
    alps::input::hybridization_parameters const p = cthyb(cthyb_required + "U_MATRIX = \"u.dat\"\n");
    EXPECT_EQ(p.u_matrix, "u.dat");
    EXPECT_FALSE(p.u_ignored);
    EXPECT_TRUE(cthyb(cthyb_required + "U_MATRIX = \"u.dat\"; U = 3\n").u_ignored);
    EXPECT_TRUE(cthyb(cthyb_required + "U = 3; MU_VECTOR = \"mu.dat\"; MU = 1\n").mu_ignored);
    EXPECT_FALSE(cthyb(cthyb_required + "U = 3; MU = 1\n").mu_ignored);
}

TEST(HybridizationInput, MeasurementsNeedTheirSizes) {
    std::string const base = cthyb_required + "U = 3\n";
    EXPECT_THROW(cthyb(base + "N_TAU = 0\n"), std::invalid_argument);
    EXPECT_THROW(cthyb(base + "MEASURE_freq = 1\n"), std::invalid_argument);
    EXPECT_THROW(cthyb(base + "MEASURE_legendre = 1; N_MATSUBARA = 10\n"), std::invalid_argument);
    EXPECT_THROW(cthyb(base + "MEASURE_nnt = 1\n"), std::invalid_argument);
    EXPECT_THROW(cthyb(base + "MEASURE_nnw = 1\n"), std::invalid_argument);
    EXPECT_THROW(cthyb(base + "MEASURE_g2w = 1; N_W = 2; N_w2 = 3\n"), std::invalid_argument);
    EXPECT_THROW(cthyb(base + "COMPUTE_VERTEX = 1; MEASURE_g2w = 1; N_W = 2; N_w2 = 4\n"), std::invalid_argument);
    EXPECT_THROW(cthyb(base + "COMPUTE_VERTEX = 1; MEASURE_freq = 1; N_MATSUBARA = 2; MEASURE_g2w = 1; N_W = 2; N_w2 = 4\n"),
                 std::invalid_argument);
    EXPECT_NO_THROW(cthyb(base + "COMPUTE_VERTEX = 1; MEASURE_freq = 1; N_MATSUBARA = 3; MEASURE_g2w = 1; N_W = 2; N_w2 = 4\n"));
}

TEST(HybridizationInput, NegativeCountsAreRejected) {
    std::string const sizes = " BETA = 4; DELTA = \"delta.dat\"; U = 3\n";
    std::string const counts[] = {"SWEEPS", "THERMALIZATION", "N_MEAS", "N_TAU", "N_ORBITALS"};
    for (std::string const & negative : counts) {
        std::string supplied = "MAX_TIME = 60;" + sizes;
        for (std::string const & key : counts)
            supplied += key + (key == negative ? " = -1\n" : " = 10\n");
        EXPECT_THROW(cthyb(supplied), alps::params_ns::exception::value_mismatch) << negative;
    }
    std::string const base = cthyb_required + "U = 3\n";
    EXPECT_THROW(cthyb(base + "N_MATSUBARA = -1\n"), alps::params_ns::exception::value_mismatch);
    EXPECT_THROW(cthyb(base + "NUM_BINS = -1\n"), alps::params_ns::exception::value_mismatch);
}
}
