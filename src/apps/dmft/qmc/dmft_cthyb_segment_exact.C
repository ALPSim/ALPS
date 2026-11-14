// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// CT-HYB-SEGMENT against exact diagonalization. The impurity, with
// interaction U, couples to two bath levels per spin, at energies -1 and +1
// with hybridization V, so that the whole system has 64 states:
//
// - G(tau) and the density at mu = U, i.e. MU = U/2 relative to half filling,
//   agree with the exact ones;
// - the same model given as the DMFT loop gives it (interaction matrix,
//   absolute orbital energies and one Delta per flavor, all in HDF5) gives
//   the same Monte Carlo run;
// - with a retarded interaction K(tau), MU = -2K'(0) is half filling, as
//   particle-hole symmetry requires when mu is shifted by +K'(0).

#include <gtest/gtest.h>
#include <alps/testing/temporary_directory.hpp>

#include <alps/hdf5/archive.hpp>
#include <alps/hdf5/vector.hpp>
#include <alps/params.hpp>
#include <alps/solvers.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#ifdef ALPS_HAVE_MPI
#include <mpi.h>
#endif

namespace {
const double beta = 10.;
const double U = 2.;
const double V = 0.5;
const double bath[] = {-1., 1.};
const int n_bath = 2;
const int n_tau = 200;

bool is_master() {
#ifdef ALPS_HAVE_MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank == 0;
#else
  return true;
#endif
}

// The solver writes text files into the working directory.
struct working_directory {
  explicit working_directory(std::filesystem::path const & dir)
    : previous(std::filesystem::current_path()) { std::filesystem::current_path(dir); }
  ~working_directory() { std::filesystem::current_path(previous); }
  std::filesystem::path previous;
};

// Exact G(tau) of flavor 0 and the density of the impurity with chemical
// potential mu. Mode s*(1+n_bath) is the impurity of spin s, the next n_bath
// modes its bath levels.
struct exact_solution {
  explicit exact_solution(double mu) {
    const int modes = 2 * (1 + n_bath), dim = 1 << modes;
    auto annihilate = [&](int mode) {
      Eigen::MatrixXd c = Eigen::MatrixXd::Zero(dim, dim);
      for (int state = 0; state < dim; ++state)
        if (state >> mode & 1) {
          const int below = std::bitset<32>(state & ((1 << mode) - 1)).count();
          c(state ^ (1 << mode), state) = below % 2 ? -1. : 1.;
        }
      return c;
    };
    std::vector<Eigen::MatrixXd> c;
    for (int mode = 0; mode < modes; ++mode) c.push_back(annihilate(mode));
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(dim, dim);
    for (int s = 0; s < 2; ++s) {
      const int d = s * (1 + n_bath);
      H -= mu * c[d].transpose() * c[d];
      for (int k = 0; k < n_bath; ++k) {
        const Eigen::MatrixXd hop = c[d].transpose() * c[d + 1 + k];
        H += bath[k] * c[d + 1 + k].transpose() * c[d + 1 + k] + V * (hop + hop.transpose());
      }
    }
    H += U * (c[0].transpose() * c[0]) * (c[1 + n_bath].transpose() * c[1 + n_bath]);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(H);
    energies = eigen.eigenvalues().array() - eigen.eigenvalues().minCoeff();
    const Eigen::MatrixXd d = eigen.eigenvectors().transpose() * c[0] * eigen.eigenvectors();
    d2 = d.cwiseAbs2();
    partition = (-beta * energies).exp().sum();
    density = (-green(beta - 1e-12));
  }
  // G(tau) = -Tr[e^{-(beta-tau)H} d e^{-tau H} d^+]/Z
  double green(double tau) const {
    double sum = 0;
    for (int a = 0; a < energies.size(); ++a)
      for (int b = 0; b < energies.size(); ++b)
        sum += std::exp(-(beta - tau) * energies[a] - tau * energies[b]) * d2(a, b);
    return -sum / partition;
  }
  Eigen::ArrayXd energies;
  Eigen::MatrixXd d2;
  double partition, density;
};

// Delta(tau) = -V^2 sum_k e^{-tau eps_k}/(1+e^{-beta eps_k})
double delta(double tau) {
  double sum = 0;
  for (double eps : bath) sum -= V * V * std::exp(-tau * eps) / (1 + std::exp(-beta * eps));
  return sum;
}

alps::params parameters(std::string const & delta_file, double mu) {
  alps::params p;
  p.supply("BETA", std::to_string(beta));
  p.supply("U", std::to_string(U));
  p.supply("MU", std::to_string(mu));
  p.supply("FLAVORS", "2");
  p.supply("N", std::to_string(n_tau));
  p.supply("NMATSUBARA", "20");
  p.supply("MAX_TIME", "600");
  p.supply("SEED", "7");
  p.supply("cthyb.DELTA", delta_file);
  p.supply("cthyb.N_MEAS", "20");
  p.supply("cthyb.SWEEPS", "4000000");
  p.supply("cthyb.THERMALIZATION", "20000");
  p.supply("cthyb.MEASURE_freq", "0");
  return p;
}

std::string write_delta_text(std::filesystem::path const & dir) {
  const std::string file = (dir / "delta.dat").string();
  std::ofstream out(file);
  out.precision(17);
  for (int i = 0; i <= n_tau; ++i) {
    const double d = delta(beta * i / n_tau);
    out << i << ' ' << d << ' ' << d << '\n';
  }
  return file;
}

struct result {
  std::vector<double> green;
  double density[2];
};

result read_result(std::string const & file) {
  alps::hdf5::archive ar(file, "r");
  result r;
  ar >> alps::make_pvp("/G_tau/0/mean/value", r.green);
  for (int f = 0; f < 2; ++f)
    ar >> alps::make_pvp("/simulation/results/density_" + std::to_string(f) + "/mean/value", r.density[f]);
  return r;
}
}

TEST(DmftRegression, cthyb_segment_exact_diagonalization) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  const std::string output = (temporary.path() / "plain.out.h5").string();
  // MU is relative to half filling: the solver adds U/2.
  alps::solvers::cthyb_segment(parameters(write_delta_text(temporary.path()), 1.), output);
  if (!is_master()) return;

  const exact_solution exact(1. + U / 2);
  const result r = read_result(output);
  ASSERT_EQ(r.green.size(), static_cast<std::size_t>(n_tau + 1));
  for (int f = 0; f < 2; ++f)
    EXPECT_NEAR(r.density[f], exact.density, 0.006) << "flavor " << f;
  // Single tau bins are noisy; compare averages over 21 bins. Over 20 seeds
  // the largest deviations are 0.0032 for the density and 0.0071 here.
  for (int center = 20; center <= n_tau - 20; center += 20) {
    double average = 0, exact_average = 0;
    for (int t = center - 10; t <= center + 10; ++t) {
      average += r.green[t] / 21;
      exact_average += exact.green(beta * t / n_tau) / 21;
    }
    EXPECT_NEAR(average, exact_average, 0.012) << "tau " << beta * center / n_tau;
  }
  EXPECT_NEAR(r.green[0], -(1 - r.density[0]), 1e-12);
  EXPECT_NEAR(r.green[n_tau], -r.density[0], 1e-12);
}

TEST(DmftRegression, cthyb_segment_dmft_inputs) {
#ifdef ALPS_HAVE_MPI
  // Under MPI the run stops at scheduled checks, not after exactly SWEEPS.
  GTEST_SKIP() << "compares two runs update by update, which needs a serial build";
#endif
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  const std::string plain = (temporary.path() / "plain.out.h5").string();
  alps::solvers::cthyb_segment(parameters(write_delta_text(temporary.path()), 1.), plain);

  // The interaction matrix, absolute orbital energies and Delta per flavor
  // in one HDF5 file, as the DMFT loop and the tutorials provide them.
  // Every rank has its own temporary directory and writes its own inputs.
  const std::string input = (temporary.path() / "input.h5").string();
  {
    alps::hdf5::archive ar(input, "w");
    std::vector<double> d(n_tau + 1);
    for (int i = 0; i <= n_tau; ++i) d[i] = delta(beta * i / n_tau);
    ar << alps::make_pvp("/Delta_0", d) << alps::make_pvp("/Delta_1", d);
    ar << alps::make_pvp("/Umatrix", std::vector<double>{0., U, U, 0.});
    ar << alps::make_pvp("/MUvector", std::vector<double>{1. + U / 2, 1. + U / 2});
  }
  alps::params p = parameters(input, 0.);
  p.supply("cthyb.DELTA_IN_HDF5", "1");
  p.supply("U_MATRIX", input);
  p.supply("UMATRIX_IN_HDF5", "1");
  p.supply("MU_VECTOR", input);
  p.supply("MU_IN_HDF5", "1");
  const std::string from_hdf5 = (temporary.path() / "hdf5.out.h5").string();
  alps::solvers::cthyb_segment(p, from_hdf5);
  if (!is_master()) return;

  const result a = read_result(plain), b = read_result(from_hdf5);
  for (int f = 0; f < 2; ++f) EXPECT_EQ(a.density[f], b.density[f]);
  EXPECT_EQ(a.green, b.green);
}

TEST(DmftRegression, cthyb_segment_retarded_interaction_half_filling) {
  alps::testing::TemporaryDirectory temporary;
  working_directory cd(temporary.path());
  // K(tau) of a Holstein phonon of frequency w0 and coupling lambda; K'(0) = lambda^2/w0
  const double lambda = 1., w0 = 5.;
  const std::string kfile = (temporary.path() / "K.h5").string();
  {
    std::vector<double> K(n_tau + 1), Kp(n_tau + 1);
    for (int i = 0; i <= n_tau; ++i) {
      const double tau = beta * i / n_tau, s = std::sinh(w0 * beta / 2);
      K[i] = -lambda * lambda / (w0 * w0) * (std::cosh(w0 * (beta / 2 - tau)) / s - std::cosh(w0 * beta / 2) / s);
      Kp[i] = lambda * lambda / w0 * std::sinh(w0 * (beta / 2 - tau)) / s;
    }
    alps::hdf5::archive ar(kfile, "w");
    ar << alps::make_pvp("/Ret_int_K", K) << alps::make_pvp("/Ret_int_Kp", Kp);
  }
  // With mu -> mu-K'(0) instead of mu+K'(0), the density is 0.44; over 10
  // seeds the largest deviation from 1/2 is 0.011 for one flavor.
  alps::params p = parameters(write_delta_text(temporary.path()), -2 * lambda * lambda / w0);
  p.supply("cthyb.RET_INT_K", kfile);
  p.supply("cthyb.K_IN_HDF5", "1");
  const std::string output = (temporary.path() / "retarded.out.h5").string();
  alps::solvers::cthyb_segment(p, output);
  if (!is_master()) return;

  const result r = read_result(output);
  EXPECT_NEAR((r.density[0] + r.density[1]) / 2, 0.5, 0.02);
}
