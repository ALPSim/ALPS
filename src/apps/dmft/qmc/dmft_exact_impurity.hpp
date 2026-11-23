// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Shared by the CT-QMC tests that compare with exact diagonalization: an
// impurity with interaction U that couples to two bath levels per spin, at
// energies -1 and +1 with hybridization V, so that the whole system has 64
// states.

#ifndef ALPS_DMFT_EXACT_IMPURITY_HPP
#define ALPS_DMFT_EXACT_IMPURITY_HPP

#include "dmft_hypothesis.hpp"

#include <alps/config.h>  // ALPS_HAVE_MPI, for is_master()

#include <Eigen/Dense>

#include <bitset>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#ifdef ALPS_HAVE_MPI
#include <mpi.h>
#endif

namespace exact_impurity {
const double beta = 10.;
const double U = 2.;
const double V = 0.5;
const double bath[] = {-1., 1.};
const int n_bath = 2;

inline bool is_master() {
#ifdef ALPS_HAVE_MPI
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank == 0;
#else
  return true;
#endif
}

// The solvers write text files into the working directory.
struct working_directory {
  explicit working_directory(std::filesystem::path const & dir)
    : previous(std::filesystem::current_path()) { std::filesystem::current_path(dir); }
  ~working_directory() { std::filesystem::current_path(previous); }
  std::filesystem::path previous;
};

// Exact G(tau) of flavor 0, the density and <n_up n_down> of the impurity
// with chemical potential mu. Mode s*(1+n_bath) is the impurity of spin s, the next n_bath
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
    const Eigen::ArrayXd weights = (-beta * energies).exp();
    partition = weights.sum();
    density = (-green(beta - 1e-12));
    const Eigen::MatrixXd nn = (c[0].transpose() * c[0]) * (c[1 + n_bath].transpose() * c[1 + n_bath]);
    double_occupancy =
        (weights * (eigen.eigenvectors().transpose() * nn * eigen.eigenvectors()).diagonal().array()).sum() / partition;
  }
  // G(tau) = -Tr[e^{-(beta-tau)H} d e^{-tau H} d^+]/Z
  double green(double tau) const {
    double sum = 0;
    for (int a = 0; a < energies.size(); ++a)
      for (int b = 0; b < energies.size(); ++b)
        sum += std::exp(-(beta - tau) * energies[a] - tau * energies[b]) * d2(a, b);
    return -sum / partition;
  }
  // G averaged over [t - 1/2, t + 1/2] of n_tau steps, as a CT-HYB solver
  // bins it at interior point t
  double green_bin(int t, int n_tau) const {
    double sum = 0;
    for (int k = 0; k < 20; ++k) sum += green(beta * (t - 0.5 + (k + 0.5) / 20) / n_tau);
    return sum / 20;
  }
  Eigen::ArrayXd energies;
  Eigen::MatrixXd d2;
  double partition, density, double_occupancy;
};

// Delta(tau) = -V^2 sum_k e^{-tau eps_k}/(1+e^{-beta eps_k})
inline double delta(double tau) {
  double sum = 0;
  for (double eps : bath) sum -= V * V * std::exp(-tau * eps) / (1 + std::exp(-beta * eps));
  return sum;
}

// Delta on n_tau+1 points, "index Delta_up Delta_down" per line
inline std::string write_delta_text(std::filesystem::path const & dir, int n_tau) {
  const std::string file = (dir / "delta.dat").string();
  std::ofstream out(file);
  out.precision(17);
  for (int i = 0; i <= n_tau; ++i) {
    const double d = delta(beta * i / n_tau);
    out << i << ' ' << d << ' ' << d << '\n';
  }
  return file;
}
}

#endif
