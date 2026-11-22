// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// ImpuritySolver: the file and command-line interface to the external
// impurity solver. The "solvers" here are shell commands that copy a prepared
// output file, so that the interface is tested without a solver.

#include "gtest/gtest.h"
#include <alps/testing/temporary_directory.hpp>
#include "cluster_transform.h"
#include "impurity_solver.h"
#include "alps/params.hpp"

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

alps::params parameters(int nmatsubara) {
  alps::params p;
  define_dca_parameters(p);
  define_solver_parameters(p);
  p["dca.LATTICE"] = "Betts2D-4A";
  p["dca.LATTICE_LIBRARY"] = DCA_CLUSTER_LIBRARY;
  p["dca.L"] = 4;
  p["dca.SITES"] = 4;
  p["MU"] = 0.3;
  p["t"] = 1.;
  p["tprime"] = 0.;
  p["U"] = 0.;
  p["BETA"] = 5.;
  p["FLAVORS"] = 2;
  p["NMATSUBARA"] = nmatsubara;
  p["N"] = 20;
  p["dca.G0_OMEGA_INPUT"] = "none";
  p["dca.G0_OMEGA_OUTPUT"] = "none";
  p["solver.INFILE_H5GF"] = "in.h5";
  p["solver.OUTFILE_H5GF"] = "out.h5";
  return p;
}

// G0 of the noninteracting lattice, and its Fourier transform
struct green_functions {
  cluster_matsubara_kspace_gf G_omega;
  cluster_itime_kspace_gf G_tau;
};
green_functions noninteracting(int nmatsubara, int ntau = 20) {
  NormalStateClusterTransformer T(parameters(nmatsubara));
  green_functions g{T.make_matsubara_gf(), T.make_itime_gf(ntau)};
  T.compute_noninteracting_G0(g.G_omega);
  alps::gf::fourier_frequency_to_time(g.G_omega, g.G_tau);
  return g;
}

void write_output(const std::string& file, const green_functions& g) {
  alps::hdf5::archive ar(file, "w");
  g.G_omega.save(ar, "/G_omega");
  g.G_tau.save(ar, "/G_tau");
}

std::string solve_error(ImpuritySolver& solver, green_functions& g) {
  green_functions out = g;
  try {
    solver.solve(g.G_omega, g.G_tau, out.G_omega, out.G_tau);
  } catch (std::exception const& e) {
    return e.what();
  }
  return "";
}

class ImpuritySolverTest : public ::testing::Test {
protected:
  alps::testing::TemporaryDirectory temporary;
  working_directory cd{temporary.path()};
};

// A "solver" that copies ready.h5 to the output file.
const std::string copy_solver = "cp ready.h5 out.h5 #";

// One that also writes its arguments, one per line, to args.txt.
const std::string recording_solver = "sh -c 'printf \"%s\\n\" \"$@\" > args.txt; cp ready.h5 out.h5' sh";

std::vector<std::string> lines(const std::string& file) {
  std::ifstream in(file);
  std::vector<std::string> result;
  for (std::string line; std::getline(in, line);) result.push_back(line);
  return result;
}
}

TEST_F(ImpuritySolverTest, ReadsTheOutputAndRemovesTheFiles) {
  green_functions g = noninteracting(16);
  write_output("ready.h5", g);
  ImpuritySolver solver(copy_solver, "params.ini", parameters(16));
  green_functions out = noninteracting(16);
  out.G_omega.initialize();
  solver.solve(g.G_omega, g.G_tau, out.G_omega, out.G_tau);
  EXPECT_NEAR((out.G_omega - g.G_omega).norm(), 0., 1e-14);
  EXPECT_FALSE(std::filesystem::exists("in.h5"));
  EXPECT_FALSE(std::filesystem::exists("out.h5"));
}

// A solver that exits normally without writing its output: an output file
// left from an earlier run must not be read as its result (#459).
TEST_F(ImpuritySolverTest, DoesNotReadAnOutputFileLeftFromBefore) {
  green_functions g = noninteracting(16);
  write_output("out.h5", g);
  ImpuritySolver solver("true", "params.ini", parameters(16));
  const std::string error = solve_error(solver, g);
  EXPECT_NE(error.find("failed to write the output file named out.h5"), std::string::npos) << error;
}

// The solver's G has to be on the meshes of the loop: here it was computed
// with another number of Matsubara frequencies (#459).
TEST_F(ImpuritySolverTest, RefusesOutputOnOtherMeshes) {
  write_output("ready.h5", noninteracting(24));
  green_functions g = noninteracting(16);
  ImpuritySolver solver(copy_solver, "params.ini", parameters(16));
  const std::string error = solve_error(solver, g);
  EXPECT_NE(error.find("out.h5: the solver's G_omega is not on the meshes of G0"), std::string::npos) << error;
}

// The input and output file may be the same: the solver has to find its
// input there, and dca reads its output from it.
TEST_F(ImpuritySolverTest, InputAndOutputCanBeTheSameFile) {
  green_functions g = noninteracting(16);
  write_output("ready.h5", g);
  alps::params p = parameters(16);
  p["solver.INFILE_H5GF"] = "io.h5";
  p["solver.OUTFILE_H5GF"] = "io.h5";
  ImpuritySolver solver("test -f io.h5 && cp ready.h5 io.h5 #", "params.ini", p);
  EXPECT_EQ(solve_error(solver, g), "");
}

// dca only prints the solver's G_tau, so it may be on another tau grid than
// G0_tau.
TEST_F(ImpuritySolverTest, AcceptsGTauOnAnotherTauGrid) {
  write_output("ready.h5", noninteracting(16, 40));
  green_functions g = noninteracting(16);
  ImpuritySolver solver(copy_solver, "params.ini", parameters(16));
  EXPECT_EQ(solve_error(solver, g), "");
}

// The MU passed to the solver is the one of the loop, to the last digit (it
// was rounded to 6 digits), and the parameter file name reaches the solver as
// one argument, even with a space or a quote in it (#460).
TEST_F(ImpuritySolverTest, PassesMuExactlyAndTheParameterFileAsOneArgument) {
  green_functions g = noninteracting(16);
  write_output("ready.h5", g);
  ImpuritySolver solver(recording_solver, "my par'ams.ini", parameters(16));
  const double mu = -0.43804012345678912;
  solver.adjust_mu(mu);
  green_functions out = g;
  solver.solve(g.G_omega, g.G_tau, out.G_omega, out.G_tau);
  const std::vector<std::string> args = lines("args.txt");
  ASSERT_EQ(args.size(), 2u);
  EXPECT_EQ(args[0], "my par'ams.ini");
  ASSERT_EQ(args[1].substr(0, 5), "--MU=");
  EXPECT_EQ(std::stod(args[1].substr(5)), mu) << args[1];
  // with the fewest digits that do
  solver.adjust_mu(0.3);
  solver.solve(g.G_omega, g.G_tau, out.G_omega, out.G_tau);
  EXPECT_EQ(lines("args.txt").back(), "--MU=0.3");
}
