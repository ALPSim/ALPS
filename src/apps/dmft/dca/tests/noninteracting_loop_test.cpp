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

#include <Eigen/Dense>

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

template<class GF = cluster_matsubara_kspace_gf> GF load(const fs::path& file, const std::string& path) {
  GF G;
  alps::hdf5::archive ar(file.string(), "r");
  G.load(ar, path);
  return G;
}

// Run dca in the directory dir_name and return that directory.
fs::path run_dca(const cluster_case& c, const std::string& dir_name, int max_it = 10,
                 bool write_iterations = true, const std::string& extra_parameters = "",
                 bool expect_success = true) {
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
        << "solver.OUTFILE_H5GF = alps_solver_out.h5gf\n"
        << extra_parameters;
  parms.close();
  const std::string command =
      "cd \"" + dir.string() + "\" && \"" DCA_EXECUTABLE "\" params.ini > dca.log 2>&1";
  const int status = std::system(command.c_str());
  if (expect_success)
    EXPECT_EQ(status, 0) << "dca failed, see " << (dir / "dca.log");
  else
    EXPECT_NE(status, 0) << "dca did not fail, see " << (dir / "dca.log");
  return dir;
}

int count_lines_containing(const std::string& text, const std::string& part) {
  std::istringstream in(text);
  int n = 0;
  for (std::string line; std::getline(in, line);) n += line.find(part) != std::string::npos;
  return n;
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

// dca.G0_OMEGA_OUTPUT holds the G0 for the next iteration also when the loop
// stops because it converged. With a tolerance this loose the first
// iteration converges, and no G0 was written at all (#465). That G0 is the
// coarse-grained lattice G0, not the isolated cluster's G0 the loop started
// from.
TEST(LoopOptions, WritesTheLastG0WhenConverged) {
  const cluster_case c{"Betts2D-4A", "cluster.xml", 4, 2, 32, -0.15};
  const fs::path dir = run_dca(c, "ConvergedG0Output", 10, false,
                               "sc.SELFCONSISTENCY_CONVERGENCE = 1e10\ndca.G0_OMEGA_OUTPUT = G0_out.h5\n");
  if (HasFailure()) return;
  const std::string log = read_file(dir / "dca.log");
  EXPECT_NE(log.find("dca stopped: the change of G"), std::string::npos) << log;
  EXPECT_EQ(log.find("iteration 1"), std::string::npos) << log;
  EXPECT_EQ(count_lines_containing(log, "spin symmetrization is"), 1) << log;
  ASSERT_TRUE(fs::exists(dir / "G0_out.h5"));
  const cluster_matsubara_kspace_gf G0_out = load(dir / "G0_out.h5", "/G0");
  const cluster_matsubara_kspace_gf G0_start = load(dir / "G0_omega.h5", "/G0");
  EXPECT_GT((G0_out - G0_start).norm(), 1e-3);
}

// In Nambu space the spin symmetrization assumes no Zeeman field; with one,
// dca refuses to run instead of removing the field from G0 in every iteration.
TEST(LoopOptions, RefusesNambuFieldWithSpinSymmetrization) {
  const cluster_case c{"Betts2D-4A", "cluster.xml", 4, 2, 32, -0.15};
  for (const char* field : {"dca.H", "dca.H_INIT"}) {
    const fs::path dir = run_dca(c, std::string("NambuField_") + field, 10, true,
                                 std::string("sc.NAMBU = true\n") + field + " = 0.3\n", false);
    EXPECT_NE(read_file(dir / "dca.log").find("needs sc.SPIN_SYMMETRIZATION = false"), std::string::npos)
        << field;
  }
}

// ---------------------------------------------------------------------------
// Double unit cell and Nambu frameworks. Their G carries a 2x2 index: for the
// double unit cell, G(K, j, spin) is the block (K, K+Q) of the cluster
// momenta K = 2i, K+Q = 2i+1; for Nambu, G(K, a, b) with a, b = particle,
// hole. The loop is checked as for the normal state: G0 of iteration 0 is the
// isolated cluster, the self-energy vanishes, the coarse-grained lattice G is
// a fixed point after iteration 0, and the local G matches a Brillouin-zone
// integral of the same 2x2 lattice problem.

typedef cluster_nambu_matsubara_kspace_gf block_gf;
typedef Eigen::Matrix2cd block;

// Inverse lattice Green's functions, written out here independently.
// Double unit cell: k and k+Q coupled by a field f (h/2 per spin for AFM,
// the staggered chemical potential for CO).
block double_cell_inverse(complex_t z, double eps_k, double eps_kq, double f) {
  block g;
  g << z - eps_k, f, f, z - eps_kq;
  return g;
}
// Nambu: particle and hole coupled by the d-wave field eta (cos kx - cos ky).
block nambu_inverse(double omega, double eps, double eta, const std::vector<double>& k) {
  const double p = std::cos(k[0]) - std::cos(k[1]);
  block g;
  g << complex_t(mu - eps, omega), -p * eta, -p * eta, complex_t(-mu + eps, omega);
  return g;
}

struct block_case {
  cluster_case cluster;
  std::string framework;  // "AFM", "CO" or "NAMBU"
};
std::ostream& operator<<(std::ostream& os, const block_case& c) {
  return os << c.framework << "_" << c.cluster.lattice;
}

const double field_init = 0.4, field = 0.25;  // h, staggered mu or eta; initial and in the loop

class BlockLoop : public ::testing::TestWithParam<block_case> {
protected:
  fs::path run() {
    const block_case& c = GetParam();
    std::ostringstream extra;
    if (c.framework == "NAMBU")
      extra << "sc.NAMBU = true\ndca.ETA1_INIT = " << field_init << "\ndca.ETA1 = " << field << "\n";
    else if (c.framework == "AFM")
      extra << "dca.AFM_SELFCONSISTENCY = true\ndca.H_INIT = " << field_init << "\ndca.H = " << field << "\n";
    else
      extra << "dca.CO_SELFCONSISTENCY = true\ndca.MU_STAG_INIT = " << field_init
            << "\ndca.MU_STAG = " << field << "\n";
    std::ostringstream name;
    name << c;
    std::string dir = name.str();
    for (char& ch : dir)
      if (ch == '-') ch = '_';
    return run_dca(c.cluster, dir, 10, true, extra.str());
  }
  // The inverse lattice G at momentum k (and k+Q) for spin s, with field f.
  block inverse(const std::vector<double>& k, double omega, int s, double f) const {
    const block_case& c = GetParam();
    if (c.framework == "NAMBU") return nambu_inverse(omega, epsilon(k, c.cluster.tprime), f, k);
    std::vector<double> kq(k);
    for (double& x : kq) x += M_PI;
    const double coupling = c.framework == "AFM" ? (s == 0 ? 1. : -1.) * f / 2. : f;
    return double_cell_inverse(complex_t(mu, omega), epsilon(k, c.cluster.tprime), epsilon(kq, c.cluster.tprime),
                               coupling);
  }
  // The 2x2 block of G at cluster index K for spin s (double cell) or of the
  // Nambu G at K (s ignored).
  block element(const block_gf& G, alps::gf::matsubara_index w, int K, int s) const {
    block b;
    if (GetParam().framework == "NAMBU") {
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
          b(i, j) = G(w, alps::gf::momentum_index(K), alps::gf::index(i), alps::gf::index(j));
    } else {
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
          b(i, j) = G(w, alps::gf::momentum_index(K + i), alps::gf::index(j), alps::gf::index(s));
    }
    return b;
  }
  int stride() const { return GetParam().framework == "NAMBU" ? 1 : 2; }
};

TEST_P(BlockLoop, ConvergesToCoarseGrainedLatticeGreensFunction) {
  const block_case& c = GetParam();
  const int n_sites = c.cluster.sites, dim = c.cluster.dimension;
  const fs::path dir = run();
  if (HasFailure()) return;
  for (int it = 0; it < 3; ++it)
    ASSERT_TRUE(fs::exists(dir / ("Greens_" + std::to_string(it) + ".h5"))) << "iteration " << it;
  ASSERT_FALSE(fs::exists(dir / "Greens_3.h5")) << "the loop did not stop at the fixed point";

  // Iteration 0: the isolated cluster with the initial field.
  const block_gf G0_isolated = load<block_gf>(dir / "Greens_0.h5", "/G0_komega");
  ASSERT_EQ(G0_isolated.mesh2().extent(), n_sites);
  for (int K = 0; K < n_sites; K += stride()) {
    std::vector<double> k(dim);
    for (int d = 0; d < dim; ++d) k[d] = G0_isolated.mesh2().points()(K, d);
    if (stride() == 2) {
      // dca pairs the cluster momenta 2i and 2i+1 as K and K+Q.
      for (int d = 0; d < dim; ++d)
        ASSERT_NEAR(std::remainder(G0_isolated.mesh2().points()(K + 1, d) - k[d] - M_PI, 2 * M_PI), 0., 1e-12)
            << "cluster momenta " << K << " and " << K + 1 << " are not K and K+Q";
    }
    for (alps::gf::matsubara_index w(0); w < n_matsubara; ++w)
      for (int s = 0; s < 2; ++s) {
        const block expected = inverse(k, G0_isolated.mesh1().points()[w()], s, field_init).inverse();
        ASSERT_LT((element(G0_isolated, w, K, s) - expected).norm(), 1e-12) << "K=" << K << " w=" << w();
      }
  }

  // Without interaction, the self-energy vanishes.
  for (int it = 0; it < 3; ++it) {
    const fs::path file = dir / ("Greens_" + std::to_string(it) + ".h5");
    const block_gf sigma = load<block_gf>(file, "/Sigma_komega");
    double largest = 0.;
    for (const complex_t* x = sigma.data().data(); x != sigma.data().data() + sigma.data().num_elements(); ++x)
      largest = std::max(largest, std::abs(*x));
    ASSERT_LT(largest, 1e-10) << "iteration " << it;
  }

  // From iteration 1 on, G0 = G is the coarse-grained lattice G with the loop's
  // field, and its average over the cluster is the local lattice G.
  const block_gf G0_1 = load<block_gf>(dir / "Greens_1.h5", "/G0_komega");
  const block_gf G_2 = load<block_gf>(dir / "Greens_2.h5", "/G_komega");
  for (std::size_t i = 0; i < G_2.data().num_elements(); ++i)
    ASSERT_LT(std::abs(G_2.data().data()[i] - G0_1.data().data()[i]), 1e-12);

  // Every element of every block, against the patch average of the 2x2
  // lattice G over dca's lattice momenta, set up here with dca's ClusterHelper
  // for the same cluster and dca.L. Only that k grid is shared with dca; the
  // dispersion and fields are the test's own. This checks the off-diagonal
  // elements, whose sign the zone average below cannot see (it vanishes for
  // the d-wave F).
  {
    alps::params p;
    p["dca.LATTICE"] = c.cluster.lattice;
    p["dca.LATTICE_LIBRARY"] = std::string(DCA_CLUSTER_DIR) + "/" + c.cluster.library;
    p["dca.SITES"] = c.cluster.sites;
    p["dca.L"] = c.cluster.L;
    p["dca.AFM_SELFCONSISTENCY"] = false;
    const ClusterHelper helper(p);
    const int rows = helper.n_lattice_momenta(), columns = dim + 1;  // the momentum and its weight
    std::vector<double> lattice(std::size_t(rows) * columns);
    for (int l = 0; l < rows; ++l)
      for (int d = 0; d < columns; ++d) lattice[l * columns + d] = helper.lattice_momentum(l, d);
    for (alps::gf::matsubara_index w(0); w < n_matsubara; w += 7) {
      const double omega = G_2.mesh1().points()[w()];
      for (int K = 0; K < n_sites; K += stride()) {
        for (int s = 0; s < 2; ++s) {
          block reference = block::Zero();
          double weight = 0.;
          std::vector<double> k(dim);
          for (int l = 0; l < rows; ++l) {
            for (int d = 0; d < dim; ++d) k[d] = G_2.mesh2().points()(K, d) + lattice[l * columns + d];
            const double w_l = lattice[l * columns + dim];
            reference += w_l * inverse(k, omega, s, field).inverse();
            weight += w_l;
          }
          reference /= weight;
          ASSERT_LT((element(G_2, w, K, s) - reference).norm(), 1e-10 * reference.norm())
              << "K=" << K << " w=" << w() << " spin " << s << "\n" << element(G_2, w, K, s) << "\nvs\n" << reference;
        }
      }
    }
  }

  const int m = reference_grid[dim];
  long n = 1;
  for (int d = 0; d < dim; ++d) n *= m;
  for (alps::gf::matsubara_index w(0); w < n_matsubara; w += 7) {
    const double omega = G_2.mesh1().points()[w()];
    for (int s = 0; s < 2; ++s) {
      // Local G: the diagonal of the 2x2 lattice G averaged over the zone.
      complex_t reference[2] = {0., 0.};
      std::vector<double> k(dim);
      for (long i = 0; i < n; ++i) {
        long j = i;
        for (int d = 0; d < dim; ++d, j /= m) k[d] = 2. * M_PI * (j % m + 0.5) / m;
        const block g = inverse(k, omega, s, field).inverse();
        reference[0] += g(0, 0);
        reference[1] += g(1, 1);
      }
      complex_t local[2] = {0., 0.};
      for (int K = 0; K < n_sites; K += stride()) {
        const block b = element(G_2, w, K, s);
        local[0] += b(0, 0);
        local[1] += b(1, 1);
      }
      for (int a = 0; a < 2; ++a) {
        reference[a] /= double(n);
        local[a] /= double(n_sites / stride());  // the number of blocks
      }
      if (stride() == 2) {
        // Diagonals of all blocks together cover every cluster momentum once.
        const complex_t all = (local[0] + local[1]) / 2.;
        EXPECT_LT(std::abs(all - reference[0]), 1e-8 * std::abs(reference[0])) << "w=" << w() << " spin " << s;
      } else {
        for (int a = 0; a < 2; ++a)
          EXPECT_LT(std::abs(local[a] - reference[a]), 1e-8 * std::abs(reference[a]))
              << "w=" << w() << " component " << a << ": " << local[a] << " vs " << reference[a];
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    Frameworks, BlockLoop,
    ::testing::Values(
        // Bipartite square-lattice clusters, which contain Q = (pi, pi)
        block_case{cluster_case{"Betts2D-4A", "cluster.xml", 4, 2, 32, -0.15}, "AFM"},
        block_case{cluster_case{"Betts2D-8A", "cluster.xml", 8, 2, 24, -0.15}, "AFM"},
        block_case{cluster_case{"Betts2D-16B", "cluster.xml", 16, 2, 16, -0.15}, "AFM"},
        block_case{cluster_case{"Betts2D-4A", "cluster.xml", 4, 2, 32, -0.15}, "CO"},
        block_case{cluster_case{"Betts2D-8A", "cluster.xml", 8, 2, 24, -0.15}, "CO"},
        // Nambu: d-wave pairing on the 1-, 4-, 8- and 16-site clusters
        block_case{cluster_case{"Betts2D-1A", "cluster.xml", 1, 2, 64, -0.15}, "NAMBU"},
        block_case{cluster_case{"Betts2D-4A", "cluster.xml", 4, 2, 32, -0.15}, "NAMBU"},
        block_case{cluster_case{"Betts2D-8A", "cluster.xml", 8, 2, 24, -0.15}, "NAMBU"},
        block_case{cluster_case{"Betts2D-16B", "cluster.xml", 16, 2, 16, -0.15}, "NAMBU"}),
    [](const ::testing::TestParamInfo<block_case>& info) {
      std::string name = info.param.framework + "_" + info.param.cluster.lattice;
      for (char& ch : name)
        if (ch == '-') ch = '_';
      return name;
    });

// With the spin symmetrization off, a Nambu Zeeman field is allowed.
TEST(LoopOptions, AllowsNambuFieldWithoutSpinSymmetrization) {
  const cluster_case c{"Betts2D-4A", "cluster.xml", 4, 2, 32, -0.15};
  const fs::path dir = run_dca(c, "NambuField_unsymmetrized", 10, true,
                               "sc.NAMBU = true\ndca.H = 0.3\nsc.SPIN_SYMMETRIZATION = false\n");
  if (HasFailure()) return;
  EXPECT_NE(read_file(dir / "dca.log").find("dca stopped: the change of G"), std::string::npos);
}

// A complex d-wave pairing field is a global phase away from a real one, and
// the G0(tau) handed to the solver is real; dca refuses dca.ETA2 != 0.
TEST(LoopOptions, RefusesNambuComplexPairingField) {
  const cluster_case c{"Betts2D-4A", "cluster.xml", 4, 2, 32, -0.15};
  for (const char* field : {"dca.ETA2", "dca.ETA2_INIT"}) {
    const fs::path dir = run_dca(c, std::string("NambuEta2_") + field, 10, true,
                                 std::string("sc.NAMBU = true\n") + field + " = 0.2\n", false);
    EXPECT_NE(read_file(dir / "dca.log").find("set dca.ETA2 (and dca.ETA2_INIT) to 0 and dca.ETA1 to |eta|"), std::string::npos)
        << field;
  }
}

} // namespace
