// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Input files of the wrong size or temperature are refused (#388, cases 6
// and 7).
//
// - G0OMEGA_INPUT and G0TAU_INPUT: green_function::read read a fixed number
//   of rows. A file with fewer frequencies left the rest at 0, and a G0 of
//   another BETA was taken as it was. write() now puts a header line with
//   BETA and the dimensions first, which read() checks; a file without one
//   is checked for its size.
// - DELTA (hybridization solver): a text file longer than N_TAU+1 was cut
//   and taken to span [0, beta], one value short kept an uninitialized value,
//   and a per-flavor HDF5 array of another length was read past its end or
//   cut.

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include "green_function.h"
#include "hybridization/hybfun.hpp"

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>

#include <complex>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


namespace {
struct working_directory {
  explicit working_directory(std::filesystem::path const & dir)
    : previous(std::filesystem::current_path()) { std::filesystem::current_path(dir); }
  ~working_directory() { std::filesystem::current_path(previous); }
  std::filesystem::path previous;
};

void expect_contains(std::string const & text, std::string const & part) {
  EXPECT_NE(text.find(part), std::string::npos) << "expected \"" << part << "\" in: " << text;
}

template <class G>
std::string read_error(G & g, std::string const & file, double beta) {
  try {
    g.read(file.c_str(), beta);
  } catch (std::exception const & e) {
    return e.what();
  }
  return "";
}

matsubara_green_function_t matsubara(unsigned n, double beta) {
  matsubara_green_function_t G(n, 2);
  for (unsigned f = 0; f < 2; ++f)
    for (unsigned w = 0; w < n; ++w)
      G(w, f) = 1. / std::complex<double>(0.1 * f, (2 * w + 1) * M_PI / beta);
  return G;
}

// the table write() wrote before, without the header line
void write_without_header(std::string const & file, matsubara_green_function_t const & G) {
  std::ofstream out(file);
  out.precision(17);
  for (unsigned w = 0; w < G.nfreq(); ++w)
    out << w << " " << G(w, 0) << " " << G(w, 1) << "\n";
}

class DmftInputFiles : public ::testing::Test {
protected:
  alps::testing::TemporaryDirectory temporary;
  working_directory cd{temporary.path()};
};
}

TEST_F(DmftInputFiles, g0_file_round_trip_and_header) {
  const matsubara_green_function_t G = matsubara(40, 10.);
  G.write("G0omega_output", 10.);
  {
    std::ifstream in("G0omega_output");
    std::string header;
    std::getline(in, header);
    EXPECT_EQ(header, "# BETA = 10 NTIME = 40 NSITE = 1 NFLAVOR = 2");
  }
  matsubara_green_function_t read(40, 2);
  read.read("G0omega_output", 10.);
  for (unsigned f = 0; f < 2; ++f)
    for (unsigned w = 0; w < 40; ++w)
      EXPECT_EQ(read(w, f), G(w, f)) << w << " " << f;

  itime_green_function_t G_tau(21, 2), read_tau(21, 2);
  for (unsigned f = 0; f < 2; ++f)
    for (unsigned t = 0; t <= 20; ++t)
      G_tau(t, f) = -0.5 + 0.01 * t - 0.1 * f;
  G_tau.write("G0tau_output", 4.);
  read_tau.read("G0tau_output", 4.);
  for (unsigned f = 0; f < 2; ++f)
    for (unsigned t = 0; t <= 20; ++t)
      EXPECT_EQ(read_tau(t, f), G_tau(t, f)) << t << " " << f;
}

// The issue's case: a G0 of BETA = 5 with 20 frequencies, read by a run with
// BETA = 10 and NMATSUBARA = 40.
TEST_F(DmftInputFiles, g0_file_of_another_run) {
  matsubara(20, 5.).write("g0_20.dat", 5.);
  matsubara_green_function_t G(40, 2);
  expect_contains(read_error(G, "g0_20.dat", 10.),
                  "Green's function file g0_20.dat: was written for BETA = 5, this run has BETA = 10.");
  expect_contains(read_error(G, "g0_20.dat", 5.),
                  "has NTIME = 20, NSITE = 1, NFLAVOR = 2; this run needs NTIME = 40 (N+1 time points, or "
                  "NMATSUBARA frequencies), NSITE = 1, NFLAVOR = 2.");
  matsubara_green_function_t one_flavor(20, 1);
  expect_contains(read_error(one_flavor, "g0_20.dat", 5.), "NFLAVOR = 1.");
}

// Files written before the header: only their size can be checked.
TEST_F(DmftInputFiles, g0_file_without_header) {
  const matsubara_green_function_t G = matsubara(20, 5.);
  write_without_header("old.dat", G);
  matsubara_green_function_t read(20, 2);
  read.read("old.dat", 5.);
  for (unsigned w = 0; w < 20; ++w)
    EXPECT_EQ(read(w, 1), G(w, 1)) << w;

  matsubara_green_function_t more(40, 2);
  expect_contains(read_error(more, "old.dat", 5.),
                  "Green's function file old.dat: has 60 entries; it needs 40 rows of 3 (a label and 2 values), 120.");
  matsubara_green_function_t fewer(10, 2);
  expect_contains(read_error(fewer, "old.dat", 5.), "has 60 entries; it needs 10 rows of 3");

  std::ofstream("bad.dat") << "0 (1,0) x\n";
  matsubara_green_function_t one(1, 2);
  expect_contains(read_error(one, "bad.dat", 5.), "entry 3 (\"x\") is not a number.");
  expect_contains(read_error(one, "missing.dat", 5.), "Green's function file missing.dat: cannot open the file.");
}

namespace {
alps::input::hybridization_parameters delta_input(std::string const & file, bool in_hdf5) {
  alps::input::hybridization_parameters input{};
  input.delta = file;
  input.delta_in_hdf5 = in_hdf5;
  input.dmft_framework = false;
  input.n_tau = 10;
  input.n_orbitals = 2;
  input.beta = 1.;
  return input;
}

std::string delta_error(alps::input::hybridization_parameters const & input) {
  try {
    hybfun Delta(input);
  } catch (std::exception const & e) {
    return e.what();
  }
  return "";
}

void write_delta(std::string const & file, int rows, std::string const & extra = "") {
  std::ofstream out(file);
  for (int i = 0; i < rows; ++i) out << i << " -0.5 -0.25\n";
  out << extra;
}
}

TEST_F(DmftInputFiles, delta_text_file) {
  write_delta("delta.dat", 11);
  const hybfun Delta(delta_input("delta.dat", false));
  for (unsigned i = 0; i <= 10; ++i) {
    EXPECT_EQ(Delta(i, 0), -0.5);
    EXPECT_EQ(Delta(i, 1), -0.25);
  }
  // the issue's case: longer than N_TAU+1
  write_delta("delta.dat", 21);
  expect_contains(delta_error(delta_input("delta.dat", false)),
                  "hybridization file delta.dat: has 63 entries; it needs 11 rows of 3 (a label and 2 values), 33.");
  // one value short
  write_delta("delta.dat", 10, "10 -0.5\n");
  expect_contains(delta_error(delta_input("delta.dat", false)), "has 32 entries");
  write_delta("delta.dat", 10, "10 -0.5 nope\n");
  expect_contains(delta_error(delta_input("delta.dat", false)), "entry 33 (\"nope\") is not a number.");
}

TEST_F(DmftInputFiles, delta_hdf5_per_flavor) {
  {
    alps::hdf5::archive ar("delta.h5", "w");
    ar << alps::make_pvp("/Delta_0", std::vector<double>(11, -0.5));
    ar << alps::make_pvp("/Delta_1", std::vector<double>(7, -0.5));
  }
  expect_contains(delta_error(delta_input("delta.h5", true)),
                  "hybridization file delta.h5: /Delta_1 has 7 points, expected N_TAU+1 = 11");
  {
    alps::hdf5::archive ar("delta.h5", "w");
    ar << alps::make_pvp("/Delta_0", std::vector<double>(21, -0.5));
    ar << alps::make_pvp("/Delta_1", std::vector<double>(21, -0.5));
  }
  expect_contains(delta_error(delta_input("delta.h5", true)), "/Delta_0 has 21 points");
}
