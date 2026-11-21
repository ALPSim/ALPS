// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The DMFT loops say why they ended (#388, case 9). The Matsubara loop ended
// at MAX_IT without a message, and a NaN in a Green's function ended every
// loop as converged, since the NaN differences compared false and were
// skipped. A loop that reaches MAX_IT still ends normally, with a message; a
// Green's function that is not finite is an error.
//
// The loops run here with a solver whose G moves by 0.1 * 2^-k in iteration
// k and a Hilbert transform that returns G as the next G0, so the difference
// between iterations k and k-1 is 0.1 * 2^-k.

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include "selfconsistency.h"

#include <alps/parameter.h>

#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
const unsigned n_tau = 16;
const unsigned n_matsubara = 8;
const double nan = std::numeric_limits<double>::quiet_NaN();

// The loops write the Green's functions of every iteration into the working
// directory.
struct working_directory {
  explicit working_directory(std::filesystem::path const & dir)
    : previous(std::filesystem::current_path()) { std::filesystem::current_path(dir); }
  ~working_directory() { std::filesystem::current_path(previous); }
  std::filesystem::path previous;
};

double shift(int k) { return 0.1 * std::pow(0.5, k); }

// Returns NaN instead of G in iteration nan_iteration (0: never).
class MatsubaraSolver : public MatsubaraImpuritySolver {
public:
  explicit MatsubaraSolver(int nan_iteration = 0) : nan_iteration_(nan_iteration) {}
  result_type solve_omega(const matsubara_green_function_t& G0, const alps::Parameters&) {
    ++calls_;
    matsubara_green_function_t G(G0);
    for (unsigned f = 0; f < G.nflavor(); ++f)
      for (unsigned w = 0; w < G.nfreq(); ++w)
        G(w, f) += calls_ == nan_iteration_ ? std::complex<double>(nan, 0.) : shift(calls_);
    itime_green_function_t G_tau(n_tau + 1, G.nflavor());
    for (unsigned f = 0; f < G.nflavor(); ++f)
      for (unsigned t = 0; t <= n_tau; ++t)
        G_tau(t, f) = -0.5;
    return result_type(G, G_tau);
  }
private:
  int nan_iteration_, calls_ = 0;
};

class MatsubaraHilbert : public FrequencySpaceHilbertTransformer {
public:
  matsubara_green_function_t operator()(const matsubara_green_function_t& G, matsubara_green_function_t&,
                                        const double, const double, const double) const { return G; }
  matsubara_green_function_t initial_G0(const alps::Parameters& p) const {
    const double beta = p["BETA"];
    matsubara_green_function_t G0(n_matsubara, 2);
    for (unsigned f = 0; f < 2; ++f)
      for (unsigned w = 0; w < n_matsubara; ++w)
        G0(w, f) = 1. / std::complex<double>(0., (2 * w + 1) * M_PI / beta);
    return G0;
  }
};

class TauSolver : public ImpuritySolver {
public:
  explicit TauSolver(int nan_iteration = 0) : nan_iteration_(nan_iteration) {}
  result_type solve(const itime_green_function_t& G0, const alps::Parameters&) {
    ++calls_;
    itime_green_function_t G(G0);
    for (unsigned f = 0; f < G.nflavor(); ++f)
      for (unsigned t = 0; t < G.ntime(); ++t)
        G(t, f) += calls_ == nan_iteration_ ? nan : shift(calls_);
    return G;
  }
private:
  int nan_iteration_, calls_ = 0;
};

class TauHilbert : public HilbertTransformer {
public:
  itime_green_function_t operator()(const itime_green_function_t& G, double, double, double) const { return G; }
  itime_green_function_t initial_G0(const alps::Parameters&) const {
    itime_green_function_t G0(n_tau + 1, 2);
    for (unsigned f = 0; f < 2; ++f)
      for (unsigned t = 0; t <= n_tau; ++t)
        G0(t, f) = -0.5;
    return G0;
  }
};

alps::Parameters parameters(int max_it, double converged) {
  alps::Parameters p;
  p["BETA"] = 5.;
  p["MU"] = 0.;
  p["H"] = 0.;
  p["U"] = 0.;
  p["N"] = n_tau;
  p["NMATSUBARA"] = n_matsubara;
  p["FLAVORS"] = 2;
  p["SITES"] = 1;
  p["SYMMETRIZATION"] = 0;
  p["MAX_IT"] = max_it;
  p["CONVERGED"] = converged;
  p["BASENAME"] = "run";
  return p;
}

class DmftConvergence : public ::testing::Test {
protected:
  // dmft sets std::cout to 20 digits; the messages print their numbers in
  // their own precision.
  void SetUp() { precision = std::cout.precision(20); }
  void TearDown() { std::cout.precision(precision); }
  std::streamsize precision;
  alps::testing::TemporaryDirectory temporary;
  working_directory cd{temporary.path()};
};

std::string matsubara_loop(int max_it, double converged, int nan_iteration = 0) {
  alps::Parameters p = parameters(max_it, converged);
  MatsubaraSolver solver(nan_iteration);
  MatsubaraHilbert hilbert;
  testing::internal::CaptureStdout();
  try {
    selfconsistency_loop_omega(p, solver, hilbert);
  } catch (...) {
    testing::internal::GetCapturedStdout();
    throw;
  }
  return testing::internal::GetCapturedStdout();
}

std::string tau_loop(bool F_loop, int max_it, double converged, int nan_iteration = 0) {
  alps::Parameters p = parameters(max_it, converged);
  TauSolver solver(nan_iteration);
  TauHilbert hilbert;
  testing::internal::CaptureStdout();
  try {
    if (F_loop)
      F_selfconsistency_loop(p, solver, hilbert);
    else
      selfconsistency_loop(p, solver, hilbert);
  } catch (...) {
    testing::internal::GetCapturedStdout();
    throw;
  }
  return testing::internal::GetCapturedStdout();
}

void expect_contains(std::string const & out, std::string const & text) {
  EXPECT_NE(out.find(text), std::string::npos) << "expected \"" << text << "\" in:\n" << out;
}
}

TEST_F(DmftConvergence, matsubara_max_it) {
  // differences 0.025 and 0.0125 after iterations 2 and 3
  const std::string out = matsubara_loop(3, 0.001);
  expect_contains(out, "DMFT loop NOT converged: MAX_IT = 3 iterations reached, last max diff 0.0125 > CONVERGED = 0.001");
  expect_contains(out, "G0OMEGA_INPUT = G0omega_output");
  EXPECT_TRUE(std::filesystem::exists("G0omega_output"));
}

TEST_F(DmftConvergence, matsubara_converged) {
  // 0.1 * 2^-5 = 0.003125 > 0.003 >= 0.1 * 2^-6
  const std::string out = matsubara_loop(100, 0.003);
  expect_contains(out, "DMFT loop converged after 6 iterations: max diff 0.0015625 <= CONVERGED = 0.003");
  EXPECT_EQ(out.find("NOT converged"), std::string::npos);
}

TEST_F(DmftConvergence, matsubara_one_iteration) {
  const std::string out = matsubara_loop(1, 0.003);
  expect_contains(out, "DMFT loop stopped after 1 iteration (MAX_IT = 1) without a convergence check");
}

TEST_F(DmftConvergence, matsubara_nan) {
  try {
    matsubara_loop(10, 0.003, 2);
    FAIL() << "a NaN in G_omega ended the loop normally";
  } catch (std::runtime_error const & e) {
    expect_contains(e.what(), "G_omega from the impurity solver is not finite");
    expect_contains(e.what(), "in iteration 2");
  }
  // the files of iteration 1 are kept, and no G0 is written for the next run
  EXPECT_TRUE(std::filesystem::exists("G_omega_1"));
  EXPECT_FALSE(std::filesystem::exists("G0omega_output"));
}

TEST_F(DmftConvergence, tau_loops_max_it_and_converged) {
  for (bool F_loop : {false, true}) {
    SCOPED_TRACE(F_loop ? "F_selfconsistency_loop" : "selfconsistency_loop");
    // differences 0.05, 0.025, 0.0125 after iterations 1, 2, 3
    expect_contains(tau_loop(F_loop, 3, 0.001),
                    "DMFT loop NOT converged: MAX_IT = 3 iterations reached, last max diff 0.0125 > CONVERGED = 0.001");
    expect_contains(tau_loop(F_loop, 100, 0.003),
                    "DMFT loop converged after 6 iterations: max diff 0.0015625 <= CONVERGED = 0.003");
  }
}

TEST_F(DmftConvergence, tau_loops_nan) {
  for (bool F_loop : {false, true}) {
    SCOPED_TRACE(F_loop ? "F_selfconsistency_loop" : "selfconsistency_loop");
    try {
      tau_loop(F_loop, 10, 0.003, 2);
      FAIL() << "a NaN in G_tau ended the loop normally";
    } catch (std::runtime_error const & e) {
      expect_contains(e.what(), "G_tau from the impurity solver is not finite");
      expect_contains(e.what(), "in iteration 2");
    }
  }
}
