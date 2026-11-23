// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The retarded interaction K(tau), K'(tau) (RET_INT_K) of the hybridization
// and ct-hyb-segment solvers. Each solver read the file without checking its
// length: a shorter HDF5 array was read past its end, a longer one, or a text
// file with more lines, was cut to N+1 values and taken to span [0, beta],
// and a text file one value short stored an uninitialized value. Only K >= 0
// was checked (and K(0) = 0 by hybridization). The interpolation read one
// past the end at |tau| = beta.

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include "retarded_interaction.h"
#include "hybridization/hybretintfun.hpp"

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {
const std::size_t n_tau = 20;
const double beta = 4.;

struct working_directory {
  explicit working_directory(std::filesystem::path const & dir)
    : previous(std::filesystem::current_path()) { std::filesystem::current_path(dir); }
  ~working_directory() { std::filesystem::current_path(previous); }
  std::filesystem::path previous;
};

// K(tau) of a Holstein phonon of frequency w0 and coupling lambda
retarded_interaction holstein(std::size_t n = n_tau) {
  const double lambda = 1., w0 = 2., s = std::sinh(w0 * beta / 2);
  retarded_interaction r;
  for (std::size_t i = 0; i <= n; ++i) {
    const double tau = beta * i / n;
    r.K.push_back(-lambda * lambda / (w0 * w0) * (std::cosh(w0 * (beta / 2 - tau)) - std::cosh(w0 * beta / 2)) / s);
    r.Kp.push_back(lambda * lambda / w0 * std::sinh(w0 * (beta / 2 - tau)) / s);
  }
  return r;
}

void write_text(std::string const & file, retarded_interaction const & r, std::string const & extra = "") {
  std::ofstream out(file);
  out.precision(17);
  for (std::size_t i = 0; i < r.K.size(); ++i)
    out << beta * i / (r.K.size() - 1) << " " << r.K[i] << " " << r.Kp[i] << "\n";
  out << extra;
}

void write_hdf5(std::string const & file, std::vector<double> const & K, std::vector<double> const & Kp) {
  alps::hdf5::archive ar(file, "w");
  ar << alps::make_pvp("/Ret_int_K", K) << alps::make_pvp("/Ret_int_Kp", Kp);
}

std::string error_of(std::string const & file, bool in_hdf5) {
  try {
    read_retarded_interaction(file, in_hdf5, n_tau);
  } catch (std::exception const & e) {
    return e.what();
  }
  return "";
}

void expect_contains(std::string const & text, std::string const & part) {
  EXPECT_NE(text.find(part), std::string::npos) << "expected \"" << part << "\" in: " << text;
}

class DmftRetardedInteraction : public ::testing::Test {
protected:
  alps::testing::TemporaryDirectory temporary;
  working_directory cd{temporary.path()};
};
}

TEST_F(DmftRetardedInteraction, reads_text_and_hdf5) {
  const retarded_interaction k = holstein();
  write_text("K.dat", k);
  write_hdf5("K.h5", k.K, k.Kp);
  for (bool in_hdf5 : {false, true}) {
    const retarded_interaction r = read_retarded_interaction(in_hdf5 ? "K.h5" : "K.dat", in_hdf5, n_tau);
    ASSERT_EQ(r.K.size(), n_tau + 1);
    ASSERT_EQ(r.Kp.size(), n_tau + 1);
    for (std::size_t i = 0; i <= n_tau; ++i) {
      EXPECT_EQ(r.K[i], k.K[i]) << i;
      EXPECT_EQ(r.Kp[i], k.Kp[i]) << i;
    }
  }
}

TEST_F(DmftRetardedInteraction, refuses_a_wrong_length) {
  const retarded_interaction shorter = holstein(n_tau - 5), longer = holstein(2 * n_tau);
  write_hdf5("K.h5", shorter.K, shorter.Kp);
  expect_contains(error_of("K.h5", true), "/Ret_int_K has 16 and /Ret_int_Kp 16 values; with N = 20 both need N+1 = 21.");
  write_hdf5("K.h5", longer.K, longer.Kp);
  expect_contains(error_of("K.h5", true), "/Ret_int_K has 41");
  write_hdf5("K.h5", holstein().K, longer.Kp);
  expect_contains(error_of("K.h5", true), "/Ret_int_K has 21 and /Ret_int_Kp 41 values");

  write_text("K.dat", longer);
  expect_contains(error_of("K.dat", false), "has 123 numbers; with N = 20 it needs N+1 = 21 lines \"tau K K'\", 63 numbers.");
  write_text("K.dat", holstein(), "4.2 0 0\n");
  expect_contains(error_of("K.dat", false), "has 66 numbers");
  // one value short
  std::ostringstream text;
  write_text("K.dat", holstein());
  {
    std::ifstream in("K.dat");
    text << in.rdbuf();
  }
  std::string s = text.str();
  s = s.substr(0, s.find_last_of(' '));
  std::ofstream("K.dat") << s << "\n";
  expect_contains(error_of("K.dat", false), "has 62 numbers");
}

TEST_F(DmftRetardedInteraction, refuses_bad_values) {
  expect_contains(error_of("missing.dat", false), "RET_INT_K missing.dat: cannot open the file.");
  expect_contains(error_of("missing.h5", true), "RET_INT_K missing.h5: cannot open the file.");
  {
    alps::hdf5::archive ar("only_K.h5", "w");
    ar << alps::make_pvp("/Ret_int_K", holstein().K);
  }
  expect_contains(error_of("only_K.h5", true), "has no array /Ret_int_Kp (K_IN_HDF5 = 1).");

  std::ofstream("bad.dat") << "0 0 x\n";
  retarded_interaction k;
  expect_contains(error_of("bad.dat", false), "entry 3 (\"x\") is not a number.");

  k = holstein();
  k.K[3] = -0.1;
  write_hdf5("K.h5", k.K, k.Kp);
  expect_contains(error_of("K.h5", true), "K(tau) < 0 at index 3");

  k = holstein();
  k.K[0] = 0.01;
  write_hdf5("K.h5", k.K, k.Kp);
  expect_contains(error_of("K.h5", true), "K(0) = 0.01 and K(beta) = ");

  k = holstein();
  k.K[n_tau] = 0.01;
  write_hdf5("K.h5", k.K, k.Kp);
  expect_contains(error_of("K.h5", true), "K(beta) = 0.01; both are 0");

  k = holstein();
  k.Kp[7] = std::numeric_limits<double>::quiet_NaN();
  write_hdf5("K.h5", k.K, k.Kp);
  expect_contains(error_of("K.h5", true), "K or K' is not finite at index 7.");
}

// The interpolation at |tau| = beta read index N+1.
TEST_F(DmftRetardedInteraction, interpolation_at_beta) {
  const retarded_interaction k = holstein();
  write_hdf5("K.h5", k.K, k.Kp);
  alps::input::hybridization_parameters input{};
  input.ret_int_k = "K.h5";
  input.k_in_hdf5 = true;
  input.beta = beta;
  input.n_tau = n_tau;
  const ret_int_fun K(input);
  EXPECT_DOUBLE_EQ(K.interpolate(beta), k.K[n_tau]);
  EXPECT_DOUBLE_EQ(K.interpolate(-beta), k.K[n_tau]);
  EXPECT_DOUBLE_EQ(K.interpolate_deriv(beta), k.Kp[n_tau]);
  EXPECT_DOUBLE_EQ(K.interpolate_deriv(0.), k.Kp[0]);
  EXPECT_DOUBLE_EQ(K.interpolate(beta / 2), k.K[n_tau / 2]);
}
