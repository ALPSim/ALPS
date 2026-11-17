// Copyright (C) 2005-2026 by Emanuel Gull
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

// End-to-end test of the dca program at U=0. dca runs its self-consistency
// loop with free_solver, which returns G = G0. Starting from the isolated
// cluster G0 = 1/(i w_n + mu - eps(K)), the loop must
//  - write that G0 in iteration 0,
//  - find a zero self-energy in every iteration,
//  - replace G0 by the coarse-grained lattice Green's function after
//    iteration 0, which is then a fixed point, so that the loop converges in
//    iteration 2,
//  - and produce a cluster Green's function whose average over the cluster
//    momenta is the local lattice Green's function, computed here
//    independently on a fine momentum grid.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "alps/hdf5.hpp"

#include "cluster_transform.h"

namespace {

namespace fs = std::filesystem;
typedef std::complex<double> complex_t;

// At this high temperature the momentum sums converge quickly enough for the
// small k-grids that keep the test fast: on every cluster below, dca's local
// Green's function is within 5e-10 (relative) of the reference.
const double beta = 1.;
const double mu = 0.3;
const double t = 1.;
const int n_matsubara = 64;
// Points per direction of the reference momentum grids in 1, 2 and 3 dimensions.
const int reference_grid[4] = {0, 4096, 256, 160};

struct cluster_case {
  std::string lattice;   // name of the LATTICEGRAPH
  std::string library;   // cluster definition file in libcluster
  int sites;
  int dimension;
  int L;                 // dca.L: linear size of the lattice of clusters in the coarse-graining
  double tprime;
};

std::ostream& operator<<(std::ostream& os, const cluster_case& c) { return os << c.lattice; }

// The dispersion dca uses, written out here independently.
double epsilon(const std::vector<double>& k, double tprime) {
  double eps = 0.;
  for (double ki : k) eps -= 2. * t * std::cos(ki);
  if (k.size() == 1) eps -= 2. * tprime * std::cos(2. * k[0]);
  if (k.size() == 2) eps -= 4. * tprime * std::cos(k[0]) * std::cos(k[1]);
  if (k.size() == 3)
    eps -= 4. * tprime * (std::cos(k[0]) * std::cos(k[1]) + std::cos(k[1]) * std::cos(k[2]) +
                          std::cos(k[0]) * std::cos(k[2]));
  return eps;
}

// Local lattice Green's function at the Matsubara frequencies omegas, on an
// m^d midpoint grid of the Brillouin zone.
std::vector<complex_t> local_green_function(const std::vector<double>& omegas, int dimension,
                                            double tprime, int m) {
  std::vector<complex_t> sum(omegas.size(), 0.);
  long n = 1;
  for (int d = 0; d < dimension; ++d) n *= m;
  std::vector<double> k(dimension);
  for (long i = 0; i < n; ++i) {
    long j = i;
    for (int d = 0; d < dimension; ++d, j /= m) k[d] = 2. * M_PI * (j % m + 0.5) / m;
    const double xi = mu - epsilon(k, tprime);
    for (std::size_t w = 0; w < omegas.size(); ++w) sum[w] += 1. / complex_t(xi, omegas[w]);
  }
  for (complex_t& s : sum) s /= double(n);
  return sum;
}

cluster_matsubara_kspace_gf load(const fs::path& file, const std::string& path) {
  cluster_matsubara_kspace_gf G;
  alps::hdf5::archive ar(file.string(), "r");
  G.load(ar, path);
  return G;
}

// Run dca in the directory dir_name and return that directory.
fs::path run_dca(const cluster_case& c, const std::string& dir_name, int max_it = 10,
                 bool write_iterations = true) {
  const fs::path dir = fs::absolute(dir_name);
  fs::remove_all(dir);
  fs::create_directories(dir);
  std::ofstream parms(dir / "params.ini");
  parms << "BETA = " << beta << "\n"
        << "MU = " << mu << "\n"
        << "U = 0\n"
        << "t = " << t << "\n"
        << "tprime = " << c.tprime << "\n"
        << "N = 128\n"
        << "NMATSUBARA = " << n_matsubara << "\n"
        << "FLAVORS = 2\n"
        << "dca.LATTICE = " << c.lattice << "\n"
        << "dca.LATTICE_LIBRARY = " << DCA_CLUSTER_DIR << "/" << c.library << "\n"
        << "dca.SITES = " << c.sites << "\n"
        << "dca.L = " << c.L << "\n"
        << "dca.INITIAL_G0 = ISOLATED\n"
        << "dca.G0_OMEGA_INPUT = G0_omega.h5\n"
        << "dca.G0_OMEGA_OUTPUT = G0_omega.h5\n"
        << "sc.MAX_IT = " << max_it << "\n"
        << "sc.SELFCONSISTENCY_CONVERGENCE = 1e-8\n"
        << "sc.WRITE_ITERATIONS = " << (write_iterations ? "true" : "false") << "\n"
        << "solver.SOLVER = " << DCA_FREE_SOLVER << "\n"
        << "solver.INFILE_H5GF = alps_solver_in.h5gf\n"
        << "solver.OUTFILE_H5GF = alps_solver_out.h5gf\n";
  parms.close();
  const std::string command =
      "cd \"" + dir.string() + "\" && \"" DCA_EXECUTABLE "\" params.ini > dca.log 2>&1";
  EXPECT_EQ(std::system(command.c_str()), 0) << "dca failed, see " << (dir / "dca.log");
  return dir;
}

std::string read_file(const fs::path& file) {
  std::ifstream in(file);
  std::stringstream text;
  text << in.rdbuf();
  return text.str();
}

class NoninteractingLoop : public ::testing::TestWithParam<cluster_case> {
protected:
  fs::path run_dca() { return ::run_dca(GetParam(), GetParam().lattice); }
};

TEST_P(NoninteractingLoop, ConvergesToCoarseGrainedLatticeGreensFunction) {
  const cluster_case& c = GetParam();
  const fs::path dir = run_dca();
  if (HasFailure()) return;

  // Iterations 0 (isolated G0), 1 (coarse-grained G0) and 2 (unchanged: converged).
  for (int it = 0; it < 3; ++it)
    ASSERT_TRUE(fs::exists(dir / ("Greens_" + std::to_string(it) + ".h5"))) << "iteration " << it;
  ASSERT_FALSE(fs::exists(dir / "Greens_3.h5")) << "the loop did not stop at the fixed point";

  // Iteration 0 starts from the isolated cluster.
  const cluster_matsubara_kspace_gf G0_isolated = load(dir / "Greens_0.h5", "/G0_komega");
  ASSERT_EQ(G0_isolated.mesh1().extent(), n_matsubara);
  ASSERT_EQ(G0_isolated.mesh2().extent(), c.sites);
  for (alps::gf::momentum_index K(0); K < c.sites; ++K) {
    std::vector<double> k(c.dimension);
    for (int d = 0; d < c.dimension; ++d) k[d] = G0_isolated.mesh2().points()(K(), d);
    for (alps::gf::matsubara_index w(0); w < n_matsubara; ++w) {
      const complex_t expected = 1. / (complex_t(mu, G0_isolated.mesh1().points()[w()]) - epsilon(k, c.tprime));
      for (alps::gf::index f(0); f < 2; ++f)
        ASSERT_LT(std::abs(G0_isolated(w, K, f) - expected), 1e-12) << "K=" << K() << " w=" << w();
    }
  }

  // Without interaction, the self-energy vanishes in every iteration.
  for (int it = 0; it < 3; ++it) {
    const fs::path file = dir / ("Greens_" + std::to_string(it) + ".h5");
    const cluster_matsubara_kspace_gf sigma = load(file, "/Sigma_komega");
    const cluster_matsubara_kspace_gf G = load(file, "/G_komega");
    for (alps::gf::matsubara_index w(0); w < n_matsubara; ++w)
      for (alps::gf::momentum_index K(0); K < c.sites; ++K)
        for (alps::gf::index f(0); f < 2; ++f)
          ASSERT_LT(std::abs(sigma(w, K, f)), 1e-10 * std::abs(1. / G(w, K, f)))
              << "iteration " << it << " K=" << K() << " w=" << w();
  }

  // From iteration 1 on, G0 = G is the coarse-grained lattice Green's function.
  const cluster_matsubara_kspace_gf G0_1 = load(dir / "Greens_1.h5", "/G0_komega");
  const cluster_matsubara_kspace_gf G_2 = load(dir / "Greens_2.h5", "/G_komega");
  std::vector<double> omegas(n_matsubara);
  for (int w = 0; w < n_matsubara; ++w) omegas[w] = G_2.mesh1().points()[w];
  const std::vector<complex_t> reference =
      local_green_function(omegas, c.dimension, c.tprime, reference_grid[c.dimension]);
  for (alps::gf::matsubara_index w(0); w < n_matsubara; ++w) {
    for (alps::gf::index f(0); f < 2; ++f) {
      complex_t local = 0.;
      for (alps::gf::momentum_index K(0); K < c.sites; ++K) {
        ASSERT_LT(std::abs(G_2(w, K, f) - G0_1(w, K, f)), 1e-12) << "K=" << K() << " w=" << w();
        local += G_2(w, K, f);
      }
      local /= double(c.sites);
      const complex_t expected = reference[w()];
      EXPECT_LT(std::abs(local - expected), 1e-8 * std::abs(expected))
          << "w=" << w() << " flavor " << f() << ": " << local << " vs " << expected;
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    Clusters, NoninteractingLoop,
    ::testing::Values(
        // One-dimensional chains, including an odd one
        cluster_case{"Chain1D_1", "cluster1d.xml", 1, 1, 256, 0.2},
        cluster_case{"Chain1D_4", "cluster1d.xml", 4, 1, 64, 0.2},
        cluster_case{"Chain1D_7", "cluster1d.xml", 7, 1, 40, 0.2},
        // Square-lattice Betts clusters: bipartite and not, square and tilted
        cluster_case{"Betts2D-1A", "cluster.xml", 1, 2, 64, -0.15},
        cluster_case{"Betts2D-4A", "cluster.xml", 4, 2, 32, -0.15},
        cluster_case{"Betts2D-8A", "cluster.xml", 8, 2, 24, -0.15},
        cluster_case{"Betts2D-8C", "cluster.xml", 8, 2, 24, -0.15},
        cluster_case{"Betts2D-9A", "cluster.xml", 9, 2, 16, -0.15},
        cluster_case{"Betts2D-10A", "cluster.xml", 10, 2, 16, -0.15},
        cluster_case{"Betts2D-12A", "cluster.xml", 12, 2, 16, -0.15},
        cluster_case{"Betts2D-16B", "cluster.xml", 16, 2, 16, -0.15},
        cluster_case{"Betts2D-18A", "cluster.xml", 18, 2, 12, -0.15},
        cluster_case{"Betts2D-20A", "cluster.xml", 20, 2, 12, -0.15},
        cluster_case{"Betts2D-26A", "cluster.xml", 26, 2, 12, -0.15},
        cluster_case{"Betts2D-32A", "cluster.xml", 32, 2, 10, -0.15},
        // Simple-cubic clusters; Betts3D-1A with a fine grid, which took minutes to set up
        // before the lattice momenta were deduplicated with a hash
        cluster_case{"Betts3D-1A", "cluster3d.xml", 1, 3, 48, -0.1},
        cluster_case{"Betts3D-8A", "cluster3d.xml", 8, 3, 16, -0.1},
        cluster_case{"Kent3D-8B", "cluster3d.xml", 8, 3, 16, -0.1},
        cluster_case{"Betts3D-10A", "cluster3d.xml", 10, 3, 16, -0.1},
        cluster_case{"Betts3D-12A", "cluster3d.xml", 12, 3, 16, -0.1}),
    [](const ::testing::TestParamInfo<cluster_case>& info) {
      std::string name = info.param.lattice;
      for (char& ch : name)
        if (ch == '-') ch = '_';
      return name;
    });

// The loop needs three iterations here. Stopped earlier by sc.MAX_IT, dca says
// so and exits normally: a run for a fixed number of iterations is the usual
// way to use it, and convergence is checked afterwards.
TEST(LoopOptions, StopsAtMaxIt) {
  const cluster_case c{"Betts2D-4A", "cluster.xml", 4, 2, 32, -0.15};
  const fs::path dir = run_dca(c, "MaxIt_Betts2D-4A", 1);
  if (HasFailure()) return;
  EXPECT_TRUE(fs::exists(dir / "Greens_0.h5"));
  EXPECT_FALSE(fs::exists(dir / "Greens_1.h5"));
  EXPECT_NE(read_file(dir / "dca.log").find("dca stopped after sc.MAX_IT=1 iterations"),
            std::string::npos);
}

// sc.WRITE_ITERATIONS=false writes only the G0 to continue from; a converged
// run says that it converged.
TEST(LoopOptions, WritesNoIterationFilesWhenAskedNotTo) {
  const cluster_case c{"Betts2D-4A", "cluster.xml", 4, 2, 32, -0.15};
  const fs::path dir = run_dca(c, "NoIterationFiles_Betts2D-4A", 10, false);
  if (HasFailure()) return;
  EXPECT_TRUE(fs::exists(dir / "G0_omega.h5"));
  EXPECT_FALSE(fs::exists(dir / "Greens_0.h5"));
  EXPECT_FALSE(fs::exists(dir / "selfenergy_0"));
  EXPECT_FALSE(fs::exists(dir / "G_omega"));
  const std::string log = read_file(dir / "dca.log");
  EXPECT_NE(log.find("iteration 2: change of G"), std::string::npos);
  EXPECT_NE(log.find("dca stopped: the change of G"), std::string::npos);
  EXPECT_EQ(log.find("sc.MAX_IT=10 iterations"), std::string::npos);
}

} // namespace
