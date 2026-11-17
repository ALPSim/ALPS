// Copyright (C) 2005-2026 by Emanuel Gull
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

// sigma_interpolate_normal and sigma_interpolate_sc carry a self-energy over
// to a run with other parameters. Here the old run has NMATSUBARA = 64 and
// the new one 32 or 128 at the same BETA, so the new frequencies are old ones
// or lie above them: the interpolated self-energy must equal the old one
// where they coincide and its last value above the old range.

#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>

#include "gtest/gtest.h"

#include "alps/hdf5.hpp"
#include "alps/params.hpp"

#include "cluster_transform.h"
#include "test_dump_file.h"

namespace {

namespace fs = std::filesystem;

const double beta = 5.;
const int n_old = 64;

void write_parameters(const fs::path& file, int n_matsubara, const std::string& dump) {
  std::ofstream f(file);
  f << "BETA = " << beta << "\nMU = 0.3\nU = 2\nt = 1\nN = 100\nNMATSUBARA = " << n_matsubara << "\n"
    << "dca.LATTICE = Betts2D-4A\ndca.SITES = 4\ndca.L = 8\n"
    << "dca.LATTICE_LIBRARY = " << DCA_CLUSTER_LIBRARY << "\ndca.LATTICE_DUMP = " << dump << "\n"
    << "dca.G0_OMEGA_INPUT = G0.h5\ndca.G0_OMEGA_OUTPUT = G0.h5\n";
}

alps::params transformer_parameters(const std::string& dump) {
  alps::params p;
  p["dca.L"] = 8;
  p["dca.LATTICE"] = "Betts2D-4A";
  p["dca.SITES"] = 4;
  p["dca.LATTICE_LIBRARY"] = DCA_CLUSTER_LIBRARY;
  p["dca.LATTICE_DUMP"] = dump;
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

// Runs a tool in its own directory: interpolate from n_old to n_new.
fs::path run_tool(const std::string& tool, const std::string& name, int n_new, const std::string& extra,
                  const std::function<void(alps::hdf5::archive&, const std::string&)>& write_sigma) {
  const fs::path dir = fs::absolute(name);
  fs::remove_all(dir);
  fs::create_directories(dir);
  const std::string dump = (dir / "cluster.h5").string();
  write_parameters(dir / "old.ini", n_old, dump);
  write_parameters(dir / "new.ini", n_new, dump);
  {
    alps::hdf5::archive ar((dir / "sigma.h5").string(), "w");
    write_sigma(ar, dump);
  }
  const std::string command = "cd \"" + dir.string() + "\" && \"" + tool +
                              "\" --pfile_old=old.ini --pfile_new=new.ini --sigma=sigma.h5 --g0=out.h5 " +
                              extra + " > tool.log 2>&1";
  EXPECT_EQ(std::system(command.c_str()), 0) << "see " << (dir / "tool.log");
  return dir;
}

// The value the interpolation must give at new frequency w.
int source_frequency(int w) { return w < n_old ? w : n_old - 1; }

class InterpolationTools : public ::testing::TestWithParam<int> {};

TEST_P(InterpolationTools, NormalState) {
  const int n_new = GetParam();
  const fs::path dir = run_tool(DCA_INTERPOLATE_NORMAL, "normal_" + std::to_string(n_new), n_new, "",
                                [](alps::hdf5::archive& ar, const std::string& dump) {
    NormalStateClusterTransformer T(transformer_parameters(dump));
    cluster_matsubara_kspace_gf sigma(alps::gf::omega_k_sigma_gf(
        alps::gf::matsubara_positive_mesh(beta, n_old), T.momentum_mesh(), alps::gf::index_mesh(2)));
    for (alps::gf::matsubara_index w(0); w < n_old; ++w)
      for (alps::gf::momentum_index k(0); k < 4; ++k)
        for (alps::gf::index f(0); f < 2; ++f) sigma(w, k, f) = sigma_value(w(), k(), f(), 0);
    sigma.save(ar, "/Sigma");
  });
  if (HasFailure()) return;
  cluster_matsubara_kspace_gf result;
  alps::hdf5::archive ar((dir / "out.h5").string(), "r");
  result.load(ar, "/Sigma");
  ASSERT_EQ(result.mesh1().extent(), n_new);
  for (alps::gf::matsubara_index w(0); w < n_new; ++w)
    for (alps::gf::momentum_index k(0); k < 4; ++k)
      for (alps::gf::index f(0); f < 2; ++f)
        ASSERT_LT(std::abs(result(w, k, f) - sigma_value(source_frequency(w()), k(), f(), 0)), 1e-12)
            << "w=" << w() << " k=" << k() << " f=" << f();
}

TEST_P(InterpolationTools, Nambu) {
  const int n_new = GetParam();
  const double factor = 2.;  // mult_anomalous
  const fs::path dir = run_tool(DCA_INTERPOLATE_SC, "nambu_" + std::to_string(n_new), n_new,
                                "--mult_anomalous=2", [](alps::hdf5::archive& ar, const std::string& dump) {
    NambuClusterTransformer T(transformer_parameters(dump));
    cluster_nambu_matsubara_kspace_gf sigma(alps::gf::omega_k_sigma1_sigma2_gf(
        alps::gf::matsubara_positive_mesh(beta, n_old), T.momentum_mesh(), alps::gf::index_mesh(2),
        alps::gf::index_mesh(2)));
    for (alps::gf::matsubara_index w(0); w < n_old; ++w)
      for (alps::gf::momentum_index k(0); k < 4; ++k)
        for (alps::gf::index a(0); a < 2; ++a)
          for (alps::gf::index b(0); b < 2; ++b) sigma(w, k, a, b) = sigma_value(w(), k(), a(), b());
    sigma.save(ar, "/Sigma");
  });
  if (HasFailure()) return;
  cluster_nambu_matsubara_kspace_gf result;
  alps::hdf5::archive ar((dir / "out.h5").string(), "r");
  result.load(ar, "/Sigma");
  ASSERT_EQ(result.mesh1().extent(), n_new);
  for (alps::gf::matsubara_index w(0); w < n_new; ++w)
    for (alps::gf::momentum_index k(0); k < 4; ++k)
      for (alps::gf::index a(0); a < 2; ++a)
        for (alps::gf::index b(0); b < 2; ++b) {
          const std::complex<double> expected =
              (a() == b() ? 1. : factor) * sigma_value(source_frequency(w()), k(), a(), b());
          ASSERT_LT(std::abs(result(w, k, a, b) - expected), 1e-12)
              << "w=" << w() << " k=" << k() << " (" << a() << "," << b() << ")";
        }
}

// Fewer and more frequencies than the old run.
INSTANTIATE_TEST_SUITE_P(NewMatsubara, InterpolationTools, ::testing::Values(32, 128));

}  // namespace
