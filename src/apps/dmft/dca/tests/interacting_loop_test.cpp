// Copyright (C) 2005-2026 by Emanuel Gull
// Copyright (C) 2026 by Tobias Wolf
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

// End-to-end test of the dca loop with a self-energy. dca runs with
// constant_sigma_solver, which returns G = [G0^-1 - Sigma]^-1 for a
// K-independent Sigma(iw) = S0 + S1/(iw) chosen here. In every framework
//  - the self-energy dca extracts with Dyson's equation is that Sigma, in
//    every iteration and in every element of the 2x2 blocks;
//  - the loop reaches the fixed point after iteration 0, as at U=0;
//  - the local Green's function is that of the lattice with this Sigma,
//    integrated here independently over the Brillouin zone.
// The solver also checks that the G0_tau it receives is the Fourier transform
// of G0. With the spin symmetrization on and a Sigma that differs between the
// spins, dca must hand the spin average of the solver's G to Dyson's equation.

#include <algorithm>
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
typedef Eigen::Matrix2cd block;
typedef cluster_nambu_matsubara_kspace_gf block_gf;

const double beta = 1., mu = 0.3, t = 1., tprime = -0.15;
const int n_matsubara = 64;
const int sites = 4, L = 32;          // Betts2D-4A
const int patch_points = 24;          // Gauss-Legendre points per direction in each patch
const double field = 0.25;            // h (AFM), staggered mu (CO) or eta (Nambu)

// The self-energy the solver returns. Normal state and double cell: s0 and
// s1 on the diagonal of each spin, s0_q on the (K, K+Q) element. Nambu: s0
// for particles, -s0 for holes and s1 on both, as particle-hole symmetry of
// Sigma_up(iw) = s0 + s1/(iw) requires, and the d-wave pairing
// s0_q[0] (cos Kx - cos Ky) off the diagonal.
struct sigma_case {
  std::string framework;  // NORMAL, AFM, CO or NAMBU
  double s0[2], s0_q[2], s1;
};
std::ostream& operator<<(std::ostream& os, const sigma_case& c) { return os << c.framework; }

double epsilon(double kx, double ky) {
  return -2. * t * (std::cos(kx) + std::cos(ky)) - 4. * tprime * std::cos(kx) * std::cos(ky);
}

// Gauss-Legendre nodes and weights on [-1, 1].
std::vector<std::pair<double, double>> gauss_legendre(int n) {
  std::vector<std::pair<double, double>> rule(n);
  for (int i = 0; i < n; ++i) {
    double x = std::cos(M_PI * (i + 0.75) / (n + 0.5)), dp = 0.;
    for (int iteration = 0; iteration < 100; ++iteration) {
      double p0 = 1., p1 = x;
      for (int k = 2; k <= n; ++k) {
        const double p2 = ((2. * k - 1.) * x * p1 - (k - 1.) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      dp = n * (x * p1 - p0) / (x * x - 1.);
      const double dx = p1 / dp;
      x -= dx;
      if (std::abs(dx) < 1e-15) break;
    }
    rule[i] = {x, 2. / ((1. - x * x) * dp * dp)};
  }
  return rule;
}

std::string parameters(const sigma_case& c) {
  std::ostringstream p;
  p << "sigma_solver.S0_UP = " << c.s0[0] << "\nsigma_solver.S0_DN = " << c.s0[1]
    << "\nsigma_solver.S0_Q_UP = " << c.s0_q[0] << "\nsigma_solver.S0_Q_DN = " << c.s0_q[1]
    << "\nsigma_solver.S1_UP = " << c.s1 << "\nsigma_solver.S1_DN = " << c.s1
    << "\nsigma_solver.S0_PARTICLE = " << c.s0[0] << "\nsigma_solver.S0_HOLE = " << -c.s0[0]
    << "\nsigma_solver.S0_DWAVE = " << c.s0_q[0] << "\n";
  if (c.framework == "AFM")
    p << "dca.AFM_SELFCONSISTENCY = true\ndca.H_INIT = " << field << "\ndca.H = " << field << "\n";
  if (c.framework == "CO")
    p << "dca.CO_SELFCONSISTENCY = true\ndca.MU_STAG_INIT = " << field << "\ndca.MU_STAG = " << field << "\n";
  if (c.framework == "NAMBU")
    p << "sc.NAMBU = true\ndca.ETA1_INIT = " << field << "\ndca.ETA1 = " << field << "\n";
  return p.str();
}

fs::path run_dca(const std::string& dir_name, const std::string& extra_parameters) {
  const fs::path dir = fs::absolute(dir_name);
  fs::remove_all(dir);
  fs::create_directories(dir);
  std::ofstream parms(dir / "params.ini");
  parms << "BETA = " << beta << "\nMU = " << mu << "\nU = 2\nt = " << t << "\ntprime = " << tprime
        << "\nN = 128\nNMATSUBARA = " << n_matsubara << "\nFLAVORS = 2\n"
        << "dca.LATTICE = Betts2D-4A\ndca.LATTICE_LIBRARY = " << DCA_CLUSTER_DIR << "/cluster.xml\n"
        << "dca.SITES = " << sites << "\ndca.L = " << L << "\n"
        << "dca.INITIAL_G0 = ISOLATED\ndca.G0_OMEGA_INPUT = G0_omega.h5\ndca.G0_OMEGA_OUTPUT = G0_omega.h5\n"
        << "sc.MAX_IT = 10\nsc.SELFCONSISTENCY_CONVERGENCE = 1e-8\nsc.WRITE_ITERATIONS = true\n"
        << "solver.SOLVER = " << DCA_CONSTANT_SIGMA_SOLVER << "\n"
        << "solver.INFILE_H5GF = alps_solver_in.h5gf\nsolver.OUTFILE_H5GF = alps_solver_out.h5gf\n"
        << extra_parameters;
  parms.close();
  // The working directory is changed here rather than with "cd" in the
  // command, which would not change drives under cmd.exe.
  const fs::path cwd = fs::current_path();
  fs::current_path(dir);
  const int status = std::system("\"" DCA_EXECUTABLE "\" params.ini > dca.log 2>&1");
  fs::current_path(cwd);
  EXPECT_EQ(status, 0) << "dca failed, see " << (dir / "dca.log");
  return dir;
}

template<class GF> GF load(const fs::path& file, const std::string& path) {
  GF G;
  alps::hdf5::archive ar(file.string(), "r");
  G.load(ar, path);
  return G;
}

// The 2x2 blocks the tests compare. Normal state: G of flavor s at K on both
// diagonal elements, so that the same 2x2 algebra serves all frameworks.
// Double cell: the block of the pair (K, K+Q), K even, for spin s. Nambu: the
// block at K (s unused).
block element(const cluster_matsubara_kspace_gf& G, alps::gf::matsubara_index w, int K, int s) {
  block b = block::Zero();
  b(0, 0) = b(1, 1) = G(w, alps::gf::momentum_index(K), alps::gf::index(s));
  return b;
}
block element(const block_gf& G, alps::gf::matsubara_index w, int K, int s, bool nambu) {
  block b;
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      b(i, j) = nambu ? G(w, alps::gf::momentum_index(K), alps::gf::index(i), alps::gf::index(j))
                      : G(w, alps::gf::momentum_index(K + i), alps::gf::index(j), alps::gf::index(s));
  return b;
}

// The injected self-energy of spin s (unused in Nambu space) as a 2x2 block,
// at the cluster momentum (Kx, Ky).
block sigma(const sigma_case& c, double omega, int s, double Kx, double Ky) {
  const complex_t tail = c.s1 / complex_t(0., omega);
  block S = block::Zero();
  if (c.framework == "NAMBU") {
    S(0, 0) = c.s0[0] + tail;
    S(1, 1) = -c.s0[0] + tail;
    S(0, 1) = S(1, 0) = c.s0_q[0] * (std::cos(Kx) - std::cos(Ky));
  } else {
    S(0, 0) = S(1, 1) = c.s0[s] + tail;
    if (c.framework != "NORMAL") S(0, 1) = S(1, 0) = c.s0_q[s];
  }
  return S;
}

// The inverse lattice Green's function without Sigma at k (and k+Q).
block lattice_inverse(const sigma_case& c, double kx, double ky, double omega, int s) {
  const complex_t z(mu, omega);
  block g = block::Zero();
  if (c.framework == "NORMAL") {
    g(0, 0) = g(1, 1) = z - epsilon(kx, ky);
  } else if (c.framework == "NAMBU") {
    const double eps = epsilon(kx, ky), p = field * (std::cos(kx) - std::cos(ky));
    g << complex_t(mu - eps, omega), -p, -p, complex_t(-mu + eps, omega);
  } else {
    const double f = c.framework == "AFM" ? (s == 0 ? 1. : -1.) * field / 2. : field;
    g << z - epsilon(kx, ky), f, f, z - epsilon(kx + M_PI, ky + M_PI);
  }
  return g;
}

class InteractingLoop : public ::testing::TestWithParam<sigma_case> {};

TEST_P(InteractingLoop, ConvergesWithTheSolversSelfEnergy) {
  const sigma_case& c = GetParam();
  const bool normal = c.framework == "NORMAL", nambu = c.framework == "NAMBU";
  const int stride = normal || nambu ? 1 : 2;
  const int spins = nambu ? 1 : 2;  // a Nambu block holds both spins
  const fs::path dir = run_dca("InteractingLoop_" + c.framework, parameters(c));
  if (HasFailure()) return;
  for (int it = 0; it < 3; ++it)
    ASSERT_TRUE(fs::exists(dir / ("Greens_" + std::to_string(it) + ".h5"))) << "iteration " << it;
  ASSERT_FALSE(fs::exists(dir / "Greens_3.h5")) << "the loop did not stop at the fixed point";

  // dca's self-energy is the solver's, element by element, in every iteration.
  for (int it = 0; it < 3; ++it) {
    const fs::path file = dir / ("Greens_" + std::to_string(it) + ".h5");
    const auto S_normal = normal ? load<cluster_matsubara_kspace_gf>(file, "/Sigma_komega") : cluster_matsubara_kspace_gf();
    const auto S_block = normal ? block_gf() : load<block_gf>(file, "/Sigma_komega");
    for (alps::gf::matsubara_index w(0); w < n_matsubara; ++w) {
      const double omega = normal ? S_normal.mesh1().points()[w()] : S_block.mesh1().points()[w()];
      for (int K = 0; K < sites; K += stride)
        for (int s = 0; s < spins; ++s) {
          const block found = normal ? element(S_normal, w, K, s) : element(S_block, w, K, s, nambu);
          const auto Kv = normal ? S_normal.mesh2().points() : S_block.mesh2().points();
          const block expected = sigma(c, omega, s, Kv(K, 0), Kv(K, 1));
          ASSERT_LT((found - expected).norm(), 1e-10 * (1. + omega))
              << "iteration " << it << " K=" << K << " spin " << s << " w=" << w() << "\n"
              << found << "\nvs\n" << expected;
        }
    }
  }

  // The local G against the lattice G with Sigma, integrated over the
  // coarse-graining patches. On Betts2D-4A these are the squares of side pi
  // around the cluster momenta, and Sigma in a patch is that of its cluster
  // momentum. Gauss-Legendre in each patch converges fast, as the lattice G is
  // smooth inside a patch. The patches are taken to be those of Betts2D-4A,
  // so the cluster momenta must be its four, (0 or pi, 0 or pi).
  const fs::path last = dir / "Greens_2.h5";
  const auto G_normal = normal ? load<cluster_matsubara_kspace_gf>(last, "/G_komega") : cluster_matsubara_kspace_gf();
  const auto G_block = normal ? block_gf() : load<block_gf>(last, "/G_komega");
  const auto Kv = normal ? G_normal.mesh2().points() : G_block.mesh2().points();
  ASSERT_EQ(Kv.shape()[0], std::size_t(sites));
  {
    std::vector<int> corners;
    for (int K = 0; K < sites; ++K) {
      int corner = 0;
      for (int d = 0; d < 2; ++d) {
        const double k = std::remainder(Kv(K, d), 2. * M_PI);  // in [-pi, pi]
        ASSERT_TRUE(std::abs(k) < 1e-12 || std::abs(std::abs(k) - M_PI) < 1e-12)
            << "cluster momentum " << K << " is not one of Betts2D-4A's; the patches below assume them";
        corner += (std::abs(k) > 1. ? 1 : 0) << d;
      }
      corners.push_back(corner);
    }
    std::sort(corners.begin(), corners.end());
    ASSERT_EQ(corners, (std::vector<int>{0, 1, 2, 3})) << "the cluster momenta are not Betts2D-4A's four";
  }
  const std::vector<std::pair<double, double>> rule = gauss_legendre(patch_points);
  for (alps::gf::matsubara_index w(0); w < n_matsubara; w += 7) {
    const double omega = normal ? G_normal.mesh1().points()[w()] : G_block.mesh1().points()[w()];
    for (int s = 0; s < spins; ++s) {
      block reference = block::Zero();
      std::vector<block> patch_reference(sites, block::Zero());
      for (int K = 0; K < sites; ++K) {
        const block S = sigma(c, omega, s, Kv(K, 0), Kv(K, 1));
        for (const auto& x : rule)
          for (const auto& y : rule) {
            const double kx = Kv(K, 0) + M_PI / 2. * x.first, ky = Kv(K, 1) + M_PI / 2. * y.first;
            patch_reference[K] += x.second * y.second / 4. * block(lattice_inverse(c, kx, ky, omega, s) - S).inverse();
          }
        reference += patch_reference[K] / double(sites);
      }
      if (nambu) {
        // Every element at every K, so that the sign of the pairing counts.
        // The d-wave Sigma jumps between patches, so dca's Simpson sums over
        // each patch with dca.L = 32 are accurate to about 1e-6 here (the
        // error falls as L^-4), not to 1e-10 as for a K-independent Sigma.
        for (int K = 0; K < sites; ++K) {
          const block& expected = patch_reference[K];
          EXPECT_LT((element(G_block, w, K, s, true) - expected).norm(), 1e-5 * expected.norm())
              << "K=" << K << " w=" << w() << "\n" << element(G_block, w, K, s, true) << "\nvs\n" << expected;
        }
      }
      block local = block::Zero();
      for (int K = 0; K < sites; K += stride)
        local += normal ? element(G_normal, w, K, s) : element(G_block, w, K, s, nambu);
      local /= double(sites / stride);
      if (normal) {
        EXPECT_LT(std::abs(local(0, 0) - reference(0, 0)), 1e-8 * std::abs(reference(0, 0))) << "w=" << w() << " spin " << s;
      } else if (nambu) {
        for (int a = 0; a < 2; ++a)
          EXPECT_LT(std::abs(local(a, a) - reference(a, a)), 1e-5 * std::abs(reference(a, a)))
              << "w=" << w() << " component " << a;
      } else {
        // The diagonals of all blocks together cover every cluster momentum once.
        const complex_t all = (local(0, 0) + local(1, 1)) / 2.;
        EXPECT_LT(std::abs(all - reference(0, 0)), 1e-8 * std::abs(reference(0, 0))) << "w=" << w() << " spin " << s;
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    Frameworks, InteractingLoop,
    ::testing::Values(sigma_case{"NORMAL", {0.6, 0.6}, {0., 0.}, 0.8},
                      // A staggered moment: the (K, K+Q) element changes sign with the spin.
                      sigma_case{"AFM", {0.6, 0.6}, {0.2, -0.2}, 0.8},
                      sigma_case{"CO", {0.6, 0.6}, {0.2, 0.2}, 0.8},
                      sigma_case{"NAMBU", {0.4, 0.4}, {0.1, 0.}, 0.8}),
    [](const ::testing::TestParamInfo<sigma_case>& info) { return info.param.framework; });

// With sc.SPIN_SYMMETRIZATION (the default) and a solver whose Sigma differs
// between the spins, dca averages the solver's G over the spins before Dyson's
// equation, so its self-energy in iteration 0 is G0^-1 - <G>^-1 for both spins.
class SpinSymmetrization : public ::testing::TestWithParam<sigma_case> {};

TEST_P(SpinSymmetrization, AveragesTheSolversGreensFunctionOverSpins) {
  const sigma_case& c = GetParam();
  const bool normal = c.framework == "NORMAL";
  const int stride = normal ? 1 : 2;
  const fs::path dir = run_dca("SpinSymmetrization_" + c.framework, parameters(c));
  if (HasFailure()) return;
  const fs::path file = dir / "Greens_0.h5";
  for (alps::gf::matsubara_index w(0); w < n_matsubara; w += 9) {
    for (int K = 0; K < sites; K += stride) {
      block G0, S[2], expected = block::Zero();
      double omega;
      if (normal) {
        const auto G0_file = load<cluster_matsubara_kspace_gf>(file, "/G0_komega");
        const auto S_file = load<cluster_matsubara_kspace_gf>(file, "/Sigma_komega");
        omega = G0_file.mesh1().points()[w()];
        G0 = element(G0_file, w, K, 0);
        for (int s = 0; s < 2; ++s) S[s] = element(S_file, w, K, s);
      } else {
        const auto G0_file = load<block_gf>(file, "/G0_komega");
        const auto S_file = load<block_gf>(file, "/Sigma_komega");
        omega = G0_file.mesh1().points()[w()];
        G0 = element(G0_file, w, K, 0, false);
        for (int s = 0; s < 2; ++s) S[s] = element(S_file, w, K, s, false);
      }
      for (int s = 0; s < 2; ++s) expected += block(block(G0.inverse()) - sigma(c, omega, s, 0., 0.)).inverse() / 2.;
      expected = block(G0.inverse()) - block(expected.inverse());
      for (int s = 0; s < 2; ++s) {
        EXPECT_LT((S[s] - expected).norm(), 1e-10 * (1. + omega))
            << "K=" << K << " spin " << s << " w=" << w() << "\n" << S[s] << "\nvs\n" << expected;
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P(Frameworks, SpinSymmetrization,
                         ::testing::Values(sigma_case{"NORMAL", {0.7, 0.4}, {0., 0.}, 0.8},
                                           sigma_case{"CO", {0.7, 0.4}, {0.3, 0.1}, 0.8}),
                         [](const ::testing::TestParamInfo<sigma_case>& info) { return info.param.framework; });

}  // namespace
