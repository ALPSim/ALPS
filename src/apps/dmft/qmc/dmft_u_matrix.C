// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The interaction matrix of the DMFT loop (#388, case 5).
//
// - dmft read U_MATRIX as lines "i j U_ij", but the hybridization and
//   ct-hyb-segment solvers read the dense matrix. On a dense file it kept
//   U(0,2) and stopped, so the half-filling shift of MU was lost (n = 0.219
//   instead of 0.5 for the hybridization solver), and ct-hyb-segment, which
//   gets the matrix from dmft, ran with that interaction. A file that could
//   not be opened, and a missing U or J, were only asserted.
// - The shift that puts half filling at MU = 0 was half the sum of column 0
//   for every flavor. It is sum_j U(f,j)/2 for flavor f.

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include "U_matrix.h"
#include "externalsolver.h"
#include "green_function.h"

#include <alps/parameter.h>
#include <alps/hdf5.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct working_directory {
  explicit working_directory(std::filesystem::path const & dir)
    : previous(std::filesystem::current_path()) { std::filesystem::current_path(dir); }
  ~working_directory() { std::filesystem::current_path(previous); }
  std::filesystem::path previous;
};

// The matrix that U = 2, J = 0.25 (U' = U - 2J) builds for two bands
const std::vector<double> issue_matrix = {0, 2, 1.25, 1.5,
                                          2, 0, 1.5, 1.25,
                                          1.25, 1.5, 0, 2,
                                          1.5, 1.25, 2, 0};
// Two bands with different intraband U: the shifts differ.
const std::vector<double> two_band_matrix = {0, 4, 1, 1,
                                             4, 0, 1, 1,
                                             1, 1, 0, 2,
                                             1, 1, 2, 0};

void write(std::string const & file, std::string const & text) { std::ofstream(file) << text; }

std::string dense(std::vector<double> const & m, unsigned n) {
  std::ostringstream out;
  for (unsigned i = 0; i < n; ++i) {
    for (unsigned j = 0; j < n; ++j) out << m[i * n + j] << " ";
    out << "\n";
  }
  return out.str();
}

alps::Parameters parameters(unsigned flavors) {
  alps::Parameters p;
  p["FLAVORS"] = flavors;
  p["SITES"] = 1;
  return p;
}

std::string error_of(alps::Parameters const & p) {
  try {
    U_matrix U(p);
  } catch (std::exception const & e) {
    return e.what();
  }
  return "";
}

void expect_contains(std::string const & text, std::string const & part) {
  EXPECT_NE(text.find(part), std::string::npos) << "expected \"" << part << "\" in: " << text;
}

class DmftUMatrix : public ::testing::Test {
protected:
  alps::testing::TemporaryDirectory temporary;
  working_directory cd{temporary.path()};
};
}

TEST_F(DmftUMatrix, dense_file_is_the_matrix_of_u_and_j) {
  write("u.dat", dense(issue_matrix, 4));
  alps::Parameters p = parameters(4);
  p["U_MATRIX"] = "u.dat";
  const U_matrix from_file(p);
  alps::Parameters q = parameters(4);
  q["U"] = 2.;
  q["J"] = 0.25;
  const U_matrix from_u_and_j(q);
  for (unsigned i = 0; i < 4; ++i) {
    for (unsigned j = 0; j < 4; ++j) {
      EXPECT_EQ(from_file(i, j), issue_matrix[4 * i + j]) << i << " " << j;
      EXPECT_EQ(from_u_and_j(i, j), issue_matrix[4 * i + j]) << i << " " << j;
    }
    EXPECT_EQ(from_file.mu_shift(i), 2.375);
    EXPECT_EQ(from_u_and_j.mu_shift(i), 2.375);
  }
  EXPECT_TRUE(from_file.uniform_mu_shift());
  EXPECT_EQ(from_file.n_nonzero(), 12);
}

TEST_F(DmftUMatrix, mu_shift_per_flavor) {
  write("u.dat", dense(two_band_matrix, 4));
  alps::Parameters p = parameters(4);
  p["U_MATRIX"] = "u.dat";
  const U_matrix U(p);
  EXPECT_EQ(U.mu_shift(0), 3.);
  EXPECT_EQ(U.mu_shift(1), 3.);
  EXPECT_EQ(U.mu_shift(2), 2.);
  EXPECT_EQ(U.mu_shift(3), 2.);
  EXPECT_FALSE(U.uniform_mu_shift());

  alps::Parameters q = parameters(2);
  q["U"] = 3.;
  const U_matrix hubbard(q);
  EXPECT_EQ(hubbard.mu_shift(0), 1.5);
  EXPECT_EQ(hubbard.mu_shift(1), 1.5);
}

TEST_F(DmftUMatrix, refused_files_and_parameters) {
  alps::Parameters p = parameters(4);
  p["U_MATRIX"] = "u.dat";

  expect_contains(error_of(p), "U_MATRIX: cannot open the file u.dat.");

  write("u.dat", "0 2 1.25 1.5\n2 0 1.5 1.25\n1.25 1.5 0 2\n1.5 1.25 2\n");
  expect_contains(error_of(p), "U_MATRIX u.dat: has 15 entries; the dense matrix for FLAVORS = 4 has 16");

  write("u.dat", dense(issue_matrix, 4) + "7\n");
  expect_contains(error_of(p), "has 17 entries");

  // the format dmft read before
  write("u.dat", "0 1 2\n1 0 2\n0 2 1.25\n2 0 1.25\n");
  expect_contains(error_of(p), "The file looks like lines \"i j U_ij\"");

  write("u.dat", "0 2 x 1.5\n");
  expect_contains(error_of(p), "entry 3 (\"x\") is not a number");

  std::vector<double> m = issue_matrix;
  m[5] = 1.;
  write("u.dat", dense(m, 4));
  expect_contains(error_of(p), "the diagonal entry U(1,1) is not 0");

  m = issue_matrix;
  m[2] = 1.;
  write("u.dat", dense(m, 4));
  expect_contains(error_of(p), "U(2,0) and U(0,2) differ");

  expect_contains(error_of(parameters(2)), "FLAVORS = 2 needs U, or the interaction matrix as U_MATRIX.");
  alps::Parameters q = parameters(4);
  q["U"] = 2.;
  expect_contains(error_of(q), "FLAVORS = 4 needs J, or the interaction matrix as U_MATRIX.");
}

namespace {
// The input file dmft writes for the solver. The executable does not exist,
// so the solver fails after the input file is written (see
// dmft_externalsolver_outfile_param).
struct solver_input {
  alps::Parameters parameters;
  std::vector<double> mu_vector, u_matrix;
};

solver_input prepare(std::string const & program, alps::Parameters p) {
  p["TMPNAME"] = "solver";
  p["SC_WRITE_DELTA"] = 1;
  p["N"] = 10;
  p["NMATSUBARA"] = 10;
  p["BETA"] = 5.;
  p["t"] = 1.;
  p["N_MEAS"] = 10;  // read by the ct-hyb-segment preparation
  p["MU"] = 0.5;
  itime_green_function_t G0(11, 1, p["FLAVORS"]);
  for (unsigned f = 0; f < G0.nflavor(); ++f)
    for (unsigned i = 0; i <= 10; ++i) G0(i, f) = -0.5;
  ExternalSolver solver(std::filesystem::current_path() / "missing" / program);
  EXPECT_THROW(solver.solve(G0, p), std::exception);
  solver_input result;
  alps::hdf5::archive ar("solver.in.h5", "r");
  ar >> alps::make_pvp("/parameters", result.parameters);
  if (ar.is_data("/MUvector")) ar >> alps::make_pvp("/MUvector", result.mu_vector);
  if (ar.is_data("/Umatrix")) ar >> alps::make_pvp("/Umatrix", result.u_matrix);
  return result;
}

alps::Parameters with_file(std::vector<double> const & m) {
  write("u.dat", dense(m, 4));
  alps::Parameters p = parameters(4);
  p["U_MATRIX"] = "u.dat";
  return p;
}
}

TEST_F(DmftUMatrix, hybridization_gets_mu_per_flavor) {
  // uniform shifts: one MU, as without U_MATRIX
  solver_input in = prepare("hybridization", with_file(issue_matrix));
  EXPECT_EQ(static_cast<double>(in.parameters["MU"]), 0.5 + 2.375);
  EXPECT_FALSE(in.parameters.defined("MU_VECTOR"));

  alps::Parameters u_and_j = parameters(4);
  u_and_j["U"] = 2.;
  u_and_j["J"] = 0.25;
  EXPECT_EQ(static_cast<double>(prepare("hybridization", u_and_j).parameters["MU"]), 0.5 + 2.375);

  // shifts that differ: a MU per flavor
  in = prepare("hybridization", with_file(two_band_matrix));
  EXPECT_FALSE(in.parameters.defined("MU"));
  EXPECT_EQ(in.mu_vector, (std::vector<double>{3.5, 3.5, 2.5, 2.5}));
}

TEST_F(DmftUMatrix, cthyb_segment_gets_the_matrix_and_mu_per_flavor) {
  solver_input in = prepare("ct-hyb-segment", with_file(issue_matrix));
  EXPECT_EQ(in.u_matrix, issue_matrix);
  EXPECT_EQ(in.mu_vector, std::vector<double>(4, 0.5 + 2.375));

  in = prepare("ct-hyb-segment", with_file(two_band_matrix));
  EXPECT_EQ(in.u_matrix, two_band_matrix);
  EXPECT_EQ(in.mu_vector, (std::vector<double>{3.5, 3.5, 2.5, 2.5}));
}

// With a retarded interaction the solvers shift U -> U - 2K'(0) and
// mu -> mu + K'(0), so half filling is at sum_j U(f,j)/2 - FLAVORS K'(0).
TEST_F(DmftUMatrix, retarded_interaction_shifts_mu) {
  const double Kp0 = 0.3;
  std::vector<double> K(11), Kp(11);
  for (int i = 0; i <= 10; ++i) {
    const double x = i / 10.;
    K[i] = x * (1 - x);            // K(0) = K(beta) = 0
    Kp[i] = Kp0 * (1 - 2 * x);
  }
  {
    alps::hdf5::archive ar("K.h5", "w");
    ar << alps::make_pvp("/Ret_int_K", K) << alps::make_pvp("/Ret_int_Kp", Kp);
  }
  // ct-hyb-segment, as the solver option
  alps::Parameters p = with_file(issue_matrix);
  p["cthyb.RET_INT_K"] = "K.h5";
  p["cthyb.K_IN_HDF5"] = 1;
  EXPECT_EQ(prepare("ct-hyb-segment", p).mu_vector, std::vector<double>(4, 0.5 + 2.375 - 4 * Kp0));
  // hybridization, unprefixed, two flavors: U/2 - 2K'(0)
  alps::Parameters q = parameters(2);
  q["U"] = 3.;
  q["RET_INT_K"] = "K.h5";
  q["K_IN_HDF5"] = 1;
  EXPECT_DOUBLE_EQ(static_cast<double>(prepare("hybridization", q).parameters["MU"]), 0.5 + 1.5 - 2 * Kp0);
}
