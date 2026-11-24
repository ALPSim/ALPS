// Copyright (C) 2005-2026 by Emanuel Gull
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

// sigma_interpolate_normal and sigma_interpolate_sc carry a self-energy over
// to a run with other parameters. Here the old run has NMATSUBARA = 64 and
// the new one 32 or 128 at the same BETA, so the new frequencies are old ones
// or lie above them. The test self-energy is c + d/(i w): the tools continue
// the real part with its last value and the imaginary part as 1/w above the
// old range (#467), which reproduces it exactly at every new frequency. A self-energy
// whose beta is not the old BETA must be refused.

#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>

#include "gtest/gtest.h"

#include "alps/hdf5.hpp"
#include "alps/params.hpp"

#include "cluster_transform.h"

namespace {

namespace fs = std::filesystem;

const double beta = 5.;
const int n_old = 64;

void write_parameters(const fs::path& file, int n_matsubara, const std::string& extra = "") {
  std::ofstream f(file);
  f << extra;
  f << "BETA = " << beta << "\nMU = 0.3\nU = 2\nt = 1\nN = 100\nNMATSUBARA = " << n_matsubara << "\n"
    << "dca.LATTICE = Betts2D-4A\ndca.SITES = 4\ndca.L = 8\n"
    << "dca.LATTICE_LIBRARY = " << DCA_CLUSTER_LIBRARY << "\n"
    << "dca.G0_OMEGA_INPUT = G0.h5\ndca.G0_OMEGA_OUTPUT = G0.h5\n";
}

alps::params transformer_parameters() {
  alps::params p;
  p["dca.L"] = 8;
  p["dca.LATTICE"] = "Betts2D-4A";
  p["dca.SITES"] = 4;
  p["dca.LATTICE_LIBRARY"] = DCA_CLUSTER_LIBRARY;
  p["MU"] = 0.3;
  p["t"] = 1.;
  p["tprime"] = 0.;
  p["U"] = 2.;
  p["V"] = 0.;
  p["BETA"] = beta;
  p["FLAVORS"] = 2;
  p["NMATSUBARA"] = n_old;
  p["N"] = 100;
  p["dca.G0_OMEGA_INPUT"] = "/tmp/none";
  p["dca.G0_OMEGA_OUTPUT"] = "/tmp/none";
  p["dca.INITIAL_G0"] = "FREE";
  p["dca.AFM_SELFCONSISTENCY"] = false;
  p["dca.CO_SELFCONSISTENCY"] = false;
  p["dca.H_INIT"] = 0.;
  p["dca.H"] = 0.;
  p["dca.ETA1_INIT"] = 0.;
  p["dca.ETA2_INIT"] = 0.;
  p["dca.ETA1"] = 0.;
  p["dca.ETA2"] = 0.;
  return p;
}

// A self-energy that differs between components and frequencies.
std::complex<double> sigma_value(int w, int k, int a, int b) {
  const std::complex<double> iw(0., (2 * w + 1) * M_PI / beta);
  return (0.3 + 0.1 * k + 0.05 * a - 0.02 * b) + (1.5 + 0.1 * a * b) / iw;
}

void write_normal_sigma(alps::hdf5::archive& ar, double sigma_beta) {
  NormalStateClusterTransformer T(transformer_parameters());
  cluster_matsubara_kspace_gf sigma(alps::gf::omega_k_sigma_gf(
      alps::gf::matsubara_positive_mesh(sigma_beta, n_old), T.momentum_mesh(), alps::gf::index_mesh(2)));
  for (alps::gf::matsubara_index w(0); w < n_old; ++w)
    for (alps::gf::momentum_index k(0); k < 4; ++k)
      for (alps::gf::index f(0); f < 2; ++f) sigma(w, k, f) = sigma_value(w(), k(), f(), 0);
  sigma.save(ar, "/Sigma");
}

void write_nambu_sigma(alps::hdf5::archive& ar, double sigma_beta) {
  NambuClusterTransformer T(transformer_parameters());
  cluster_nambu_matsubara_kspace_gf sigma(alps::gf::omega_k_sigma1_sigma2_gf(
      alps::gf::matsubara_positive_mesh(sigma_beta, n_old), T.momentum_mesh(), alps::gf::index_mesh(2),
      alps::gf::index_mesh(2)));
  for (alps::gf::matsubara_index w(0); w < n_old; ++w)
    for (alps::gf::momentum_index k(0); k < 4; ++k)
      for (alps::gf::index a(0); a < 2; ++a)
        for (alps::gf::index b(0); b < 2; ++b) sigma(w, k, a, b) = sigma_value(w(), k(), a(), b());
  sigma.save(ar, "/Sigma");
}

// Runs a tool in its own directory dir, interpolating a self-energy given at
// sigma_beta from n_old to n_new frequencies, and returns its exit status.
// The working directory is changed here rather than with "cd" in the command,
// which would not change drives under cmd.exe.
int run_tool(const std::string& tool, const fs::path& dir, int n_new, const std::string& extra,
             double sigma_beta, const std::function<void(alps::hdf5::archive&, double)>& write_sigma,
             const std::string& new_parameters = "") {
  fs::remove_all(dir);
  fs::create_directories(dir);
  write_parameters(dir / "old.ini", n_old);
  write_parameters(dir / "new.ini", n_new, new_parameters);
  {
    alps::hdf5::archive ar((dir / "sigma.h5").string(), "w");
    write_sigma(ar, sigma_beta);
  }
  const std::string command = "\"" + tool +
                              "\" --pfile_old=old.ini --pfile_new=new.ini --sigma=sigma.h5 --g0=out.h5 " +
                              extra + " > tool.log 2>&1";
  const fs::path cwd = fs::current_path();
  fs::current_path(dir);
  const int status = std::system(command.c_str());
  fs::current_path(cwd);
  return status;
}

std::string read_file(const fs::path& file) {
  std::ifstream f(file);
  std::stringstream s;
  s << f.rdbuf();
  return s.str();
}

// At every new frequency w the interpolation must give the self-energy itself,
// sigma_value(w, ...), also above the old range.
class InterpolationTools : public ::testing::TestWithParam<int> {};

TEST_P(InterpolationTools, NormalState) {
  const int n_new = GetParam();
  const fs::path dir = fs::absolute("normal_" + std::to_string(n_new));
  ASSERT_EQ(run_tool(DCA_INTERPOLATE_NORMAL, dir, n_new, "", beta, write_normal_sigma), 0)
      << "see " << (dir / "tool.log");
  cluster_matsubara_kspace_gf result;
  alps::hdf5::archive ar((dir / "out.h5").string(), "r");
  result.load(ar, "/Sigma");
  ASSERT_EQ(result.mesh1().extent(), n_new);
  for (alps::gf::matsubara_index w(0); w < n_new; ++w)
    for (alps::gf::momentum_index k(0); k < 4; ++k)
      for (alps::gf::index f(0); f < 2; ++f)
        ASSERT_LT(std::abs(result(w, k, f) - sigma_value(w(), k(), f(), 0)), 1e-12)
            << "w=" << w() << " k=" << k() << " f=" << f();
}

TEST_P(InterpolationTools, Nambu) {
  const int n_new = GetParam();
  const double factor = 2.;  // mult_anomalous
  const fs::path dir = fs::absolute("nambu_" + std::to_string(n_new));
  ASSERT_EQ(run_tool(DCA_INTERPOLATE_SC, dir, n_new, "--mult_anomalous=2", beta, write_nambu_sigma), 0)
      << "see " << (dir / "tool.log");
  cluster_nambu_matsubara_kspace_gf result;
  alps::hdf5::archive ar((dir / "out.h5").string(), "r");
  result.load(ar, "/Sigma");
  ASSERT_EQ(result.mesh1().extent(), n_new);
  for (alps::gf::matsubara_index w(0); w < n_new; ++w)
    for (alps::gf::momentum_index k(0); k < 4; ++k)
      for (alps::gf::index a(0); a < 2; ++a)
        for (alps::gf::index b(0); b < 2; ++b) {
          const std::complex<double> expected =
              (a() == b() ? 1. : factor) * sigma_value(w(), k(), a(), b());
          ASSERT_LT(std::abs(result(w, k, a, b) - expected), 1e-12)
              << "w=" << w() << " k=" << k() << " (" << a() << "," << b() << ")";
        }
}

// Fewer and more frequencies than the old run.
INSTANTIATE_TEST_SUITE_P(NewMatsubara, InterpolationTools, ::testing::Values(32, 128));

// A self-energy at another beta than the old BETA lives on other frequencies.
TEST(InterpolationToolsBeta, NormalStateRefusesOtherBeta) {
  const fs::path dir = fs::absolute("normal_other_beta");
  EXPECT_NE(run_tool(DCA_INTERPOLATE_NORMAL, dir, n_old, "", 2 * beta, write_normal_sigma), 0);
  EXPECT_NE(read_file(dir / "tool.log").find("BETA of the old parameter file"), std::string::npos)
      << read_file(dir / "tool.log");
}

TEST(InterpolationToolsBeta, NambuRefusesOtherBeta) {
  const fs::path dir = fs::absolute("nambu_other_beta");
  EXPECT_NE(run_tool(DCA_INTERPOLATE_SC, dir, n_old, "", 2 * beta, write_nambu_sigma), 0);
  EXPECT_NE(read_file(dir / "tool.log").find("BETA of the old parameter file"), std::string::npos)
      << read_file(dir / "tool.log");
}

// A beta within the relative 1e-10 the tools allow is the old BETA: the check
// of the momenta and flavors (#471) must not compare it exactly.
TEST(InterpolationToolsBeta, NormalStateAcceptsBetaWithinTolerance) {
  const fs::path dir = fs::absolute("normal_close_beta");
  EXPECT_EQ(run_tool(DCA_INTERPOLATE_NORMAL, dir, 32, "", beta * (1 + 1e-12), write_normal_sigma), 0)
      << read_file(dir / "tool.log");
}

TEST(InterpolationToolsBeta, NambuAcceptsBetaWithinTolerance) {
  const fs::path dir = fs::absolute("nambu_close_beta");
  EXPECT_EQ(run_tool(DCA_INTERPOLATE_SC, dir, 32, "", beta * (1 + 1e-12), write_nambu_sigma), 0)
      << read_file(dir / "tool.log");
}

}  // namespace

// The starting G0 for a new Nambu run has the new run's initial fields in its
// data (the k-sum) and in its tail: the tail was attached with dca.ETA1 and
// dca.H instead (#469). With dca.ETA1_INIT = 0.4 and dca.ETA1 = 0, the
// off-diagonal G0 at the highest frequency, (iw)^2 G0_12 = c2_12 + O(1/w),
// has to match the tail's c2_12, the patch average of eta times the d-wave
// form factor, which is not zero for every momentum.
TEST(InterpolationTools, NambuTailHasTheInitialFields) {
  const fs::path dir = fs::absolute("nambu_tail");
  ASSERT_EQ(run_tool(DCA_INTERPOLATE_SC, dir, 1024, "", beta, write_nambu_sigma,
                     "dca.ETA1_INIT = 0.4\ndca.ETA1 = 0\nsc.NAMBU = true\n"), 0)
      << "see " << (dir / "tool.log");
  cluster_nambu_matsubara_kspace_gf G0;
  alps::hdf5::archive ar((dir / "out.h5").string(), "r");
  G0.load(ar, "/G0");
  const alps::gf::matsubara_index last(G0.mesh1().extent() - 1);
  const std::complex<double> iw(0., G0.mesh1().points()[last()]);
  double largest = 0.;
  for (alps::gf::momentum_index k(0); k < G0.mesh2().extent(); ++k) {
    const double c2 = G0.tail(2)(k, alps::gf::index(0), alps::gf::index(1));
    largest = std::max(largest, std::abs(c2));
    EXPECT_NEAR(std::real(iw * iw * G0(last, k, alps::gf::index(0), alps::gf::index(1))), c2, 1e-2) << "k=" << k();
  }
  EXPECT_GT(largest, 0.1);
}

TEST(InterpolationTools, NambuRefusesComplexPairingField) {
  const fs::path dir = fs::absolute("nambu_eta2");
  EXPECT_NE(run_tool(DCA_INTERPOLATE_SC, dir, 32, "", beta, write_nambu_sigma, "dca.ETA2_INIT = 0.1\n"), 0);
  EXPECT_NE(read_file(dir / "tool.log").find("set dca.ETA2 and dca.ETA2_INIT to 0"), std::string::npos)
      << read_file(dir / "tool.log");
}

// The self-energy has to be on the momenta of the old cluster with 2 flavors:
// load() takes the meshes from the file, and a self-energy with fewer
// momenta or flavors was read past its end (#471).
namespace {
void write_mismatched_sigma(alps::hdf5::archive& ar, const char* lattice, int sites, int flavors) {
  alps::params p = transformer_parameters();
  p["dca.LATTICE"] = lattice;
  p["dca.SITES"] = sites;
  NormalStateClusterTransformer T(p);
  cluster_matsubara_kspace_gf sigma(alps::gf::omega_k_sigma_gf(
      alps::gf::matsubara_positive_mesh(beta, n_old), T.momentum_mesh(), alps::gf::index_mesh(flavors)));
  sigma.initialize();
  sigma.save(ar, "/Sigma");
}
}

TEST(InterpolationTools, RefusesSigmaOfAnotherClusterOrFlavors) {
  struct mismatch { const char* name; const char* lattice; int sites; int flavors; };
  for (mismatch m : {mismatch{"flavors", "Betts2D-4A", 4, 1}, mismatch{"cluster", "Betts2D-8A", 8, 2}}) {
    const fs::path dir = fs::absolute(std::string("mismatch_") + m.name);
    EXPECT_NE(run_tool(DCA_INTERPOLATE_NORMAL, dir, 32, "", beta,
                       [&](alps::hdf5::archive& ar, double) { write_mismatched_sigma(ar, m.lattice, m.sites, m.flavors); }),
              0) << m.name;
    EXPECT_NE(read_file(dir / "tool.log").find("sigma.h5: the self-energy is not on the momenta of the old cluster"),
              std::string::npos) << m.name << ": " << read_file(dir / "tool.log");
  }
}

TEST(InterpolationTools, RefusesOtherFlavors) {
  const fs::path dir = fs::absolute("flavors_4");
  EXPECT_NE(run_tool(DCA_INTERPOLATE_NORMAL, dir, 32, "", beta, write_normal_sigma, "FLAVORS = 4\n"), 0);
  EXPECT_NE(read_file(dir / "tool.log").find("FLAVORS = 4: the self-energy interpolation works for FLAVORS = 2."),
            std::string::npos) << read_file(dir / "tool.log");
}
