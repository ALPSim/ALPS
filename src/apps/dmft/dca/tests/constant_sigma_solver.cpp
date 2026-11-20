// Copyright (C) 2005-2026 by Emanuel Gull
// Copyright (C) 2026 by Tobias Wolf
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

// Impurity solver for the interacting loop test. Called by dca like a real
// solver, as constant_sigma_solver <parameter file> [--MU=<mu>], it reads G0
// and G0_tau from solver.INFILE_H5GF and writes
//
//     G(K, iw) = [G0(K, iw)^-1 - Sigma(iw)]^-1,   Sigma(iw) = S0 + S1/(iw),
//
// with its tail and its Fourier transform G_tau, to solver.OUTFILE_H5GF.
// Apart from the d-wave anomalous term, Sigma does not depend on K. It is one number per flavor in the normal state
// and a 2x2 matrix in the double unit cell (per spin, in the K, K+Q basis) and
// in Nambu space (particle, hole):
//
//   sigma_solver.S0_UP, S0_DN     diagonal, normal state and double cell
//   sigma_solver.S0_Q_UP, S0_Q_DN the (K, K+Q) element, double cell
//   sigma_solver.S0_PARTICLE, S0_HOLE   diagonal, Nambu
//   sigma_solver.S0_DWAVE         the anomalous element, Nambu, times
//                                 cos Kx - cos Ky (the one K dependence)
//   sigma_solver.S1_UP, S1_DN     the 1/(iw) term on the diagonal (Nambu:
//                                 particle, hole)
//
// The solver also checks that the G0_tau it gets is the Fourier transform of
// its G0, and stops with an error if it is not.
//
// --MU, which dca passes with sc.ADJUST_MU, is accepted and ignored: the
// chemical potential enters only through G0, which dca computes with it.

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "alps/gf/fourier.hpp"
#include "alps/hdf5.hpp"
#include "alps/params.hpp"

#include "cluster_transform.h"

namespace {

namespace g = alps::gf;
typedef std::complex<double> complex_t;
typedef Eigen::MatrixXcd cmatrix;
typedef Eigen::MatrixXd rmatrix;

enum framework_t { NORMAL, DOUBLE_CELL, NAMBU };

struct sigma_t {
  double s0[2], s0_q[2], s0_nambu[2], s0_dwave, s1[2];
};

// The elements of one block of G: block (K, flavor) is 1x1 in the normal
// state; the double cell block of the pair (K, K+Q) and spin s holds
// G(K+i, j, s); the Nambu block at K holds G(K, i, j).
struct element_t { int k, i, j; };
struct block_t {
  int dim;
  int spin;     // -1 in Nambu space
  double form;  // cos Kx - cos Ky in Nambu space
  std::vector<element_t> e;  // e[i*dim+j]
};

complex_t& at(cluster_matsubara_kspace_gf& G, int w, const element_t& x) {
  return G(g::matsubara_index(w), g::momentum_index(x.k), g::index(x.i));
}
complex_t& at(cluster_nambu_matsubara_kspace_gf& G, int w, const element_t& x) {
  return G(g::matsubara_index(w), g::momentum_index(x.k), g::index(x.i), g::index(x.j));
}
double& at(cluster_matsubara_kspace_gf::tail_type& T, const element_t& x) {
  return T(g::momentum_index(x.k), g::index(x.i));
}
double& at(cluster_nambu_matsubara_kspace_gf::tail_type& T, const element_t& x) {
  return T(g::momentum_index(x.k), g::index(x.i), g::index(x.j));
}

template<class GFW, class GFT>
void solve(const std::string& infile, const std::string& outfile, framework_t framework, const sigma_t& s) {
  GFW G0;
  GFT G0_tau;
  {
    alps::hdf5::archive in(infile, "r");
    G0.load(in, "/G0");
    G0_tau.load(in, "/G0_tau");
  }
  if (G0.min_tail_order() == g::TAIL_NOT_SET || G0.min_tail_order() > 1 || G0.max_tail_order() < 3)
    throw std::runtime_error("G0 has no tail of orders 1 to 3");

  // The G0_tau dca hands over must be the transform of its G0.
  {
    GFT expected = G0_tau;
    g::fourier_frequency_to_time(G0, expected);
    double largest = 0.;
    for (std::size_t i = 0; i < expected.data().num_elements(); ++i)
      largest = std::max(largest, std::abs(expected.data().data()[i] - G0_tau.data().data()[i]));
    if (largest > 1e-10)
      throw std::runtime_error("G0_tau is not the Fourier transform of G0 (largest difference " +
                               std::to_string(largest) + ")");
  }

  const int n_freq = G0.mesh1().extent(), n_k = G0.mesh2().extent();
  std::vector<block_t> blocks;
  if (framework == NORMAL) {
    for (int k = 0; k < n_k; ++k)
      for (int f = 0; f < 2; ++f) blocks.push_back(block_t{1, f, 0., {element_t{k, f, 0}}});
  } else if (framework == DOUBLE_CELL) {
    for (int k = 0; k < n_k; k += 2)
      for (int spin = 0; spin < 2; ++spin) {
        block_t b{2, spin, 0., {}};
        for (int i = 0; i < 2; ++i)
          for (int j = 0; j < 2; ++j) b.e.push_back(element_t{k + i, j, spin});
        blocks.push_back(b);
      }
  } else {
    for (int k = 0; k < n_k; ++k) {
      const auto K = G0.mesh2().points();
      block_t b{2, -1, std::cos(K(k, 0)) - std::cos(K(k, 1)), {}};
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) b.e.push_back(element_t{k, i, j});
      blocks.push_back(b);
    }
  }

  GFW G = G0;
  typename GFW::tail_type c1(G0.tail(1)), c2(G0.tail(2)), c3(G0.tail(3));
  for (const block_t& b : blocks) {
    rmatrix S0 = rmatrix::Zero(b.dim, b.dim), S1 = rmatrix::Zero(b.dim, b.dim);
    if (framework == NORMAL) {
      S0(0, 0) = s.s0[b.spin];
      S1(0, 0) = s.s1[b.spin];
    } else if (framework == DOUBLE_CELL) {
      S0(0, 0) = S0(1, 1) = s.s0[b.spin];
      S0(0, 1) = S0(1, 0) = s.s0_q[b.spin];
      S1(0, 0) = S1(1, 1) = s.s1[b.spin];
    } else {
      for (int i = 0; i < 2; ++i) {
        S0(i, i) = s.s0_nambu[i];
        S1(i, i) = s.s1[i];
      }
      S0(0, 1) = S0(1, 0) = s.s0_dwave * b.form;
    }
    for (int w = 0; w < n_freq; ++w) {
      const complex_t iw(0., G0.mesh1().points()[w]);
      cmatrix g0(b.dim, b.dim);
      for (int i = 0; i < b.dim; ++i)
        for (int j = 0; j < b.dim; ++j) g0(i, j) = at(G0, w, b.e[i * b.dim + j]);
      const cmatrix g = (cmatrix(g0.inverse()) - S0.cast<complex_t>() - S1.cast<complex_t>() / iw).inverse();
      for (int i = 0; i < b.dim; ++i)
        for (int j = 0; j < b.dim; ++j) at(G, w, b.e[i * b.dim + j]) = g(i, j);
    }
    // With G0^-1 = iw - A - B/(iw) + ..., where A and B follow from the tail
    // of G0, G^-1 = iw - (A + S0) - (B + S1)/(iw) + ...
    rmatrix A(b.dim, b.dim), B(b.dim, b.dim);
    for (int i = 0; i < b.dim; ++i)
      for (int j = 0; j < b.dim; ++j) {
        A(i, j) = at(c2, b.e[i * b.dim + j]);
        B(i, j) = at(c3, b.e[i * b.dim + j]);
      }
    B -= A * A;
    const rmatrix C2 = A + S0, C3 = C2 * C2 + B + S1;
    for (int i = 0; i < b.dim; ++i)
      for (int j = 0; j < b.dim; ++j) {
        at(c1, b.e[i * b.dim + j]) = i == j ? 1. : 0.;
        at(c2, b.e[i * b.dim + j]) = C2(i, j);
        at(c3, b.e[i * b.dim + j]) = C3(i, j);
      }
  }
  G.set_tail(1, c1);
  G.set_tail(2, c2);
  G.set_tail(3, c3);
  GFT G_tau = G0_tau;
  g::fourier_frequency_to_time(G, G_tau);

  alps::hdf5::archive out(outfile, "w");
  G.save(out, "/G_omega");
  G_tau.save(out, "/G_tau");
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2)
      throw std::invalid_argument("usage: constant_sigma_solver <parameter file> [--MU=<mu> (ignored)]");
    alps::params p{std::string(argv[1])};
    p.define<std::string>("solver.INFILE_H5GF", "alps_solver_in.h5gf", "input file");
    p.define<std::string>("solver.OUTFILE_H5GF", "alps_solver_out.h5gf", "output file");
    p.define<bool>("sc.NAMBU", false, "Nambu space");
    for (const char* key : {"S0_UP", "S0_DN", "S0_Q_UP", "S0_Q_DN", "S0_PARTICLE", "S0_HOLE", "S0_DWAVE", "S1_UP", "S1_DN"})
      p.define<double>(std::string("sigma_solver.") + key, 0., "self-energy term");
    sigma_t s;
    s.s0[0] = p["sigma_solver.S0_UP"];
    s.s0[1] = p["sigma_solver.S0_DN"];
    s.s0_q[0] = p["sigma_solver.S0_Q_UP"];
    s.s0_q[1] = p["sigma_solver.S0_Q_DN"];
    s.s0_nambu[0] = p["sigma_solver.S0_PARTICLE"];
    s.s0_nambu[1] = p["sigma_solver.S0_HOLE"];
    s.s0_dwave = p["sigma_solver.S0_DWAVE"];
    s.s1[0] = p["sigma_solver.S1_UP"];
    s.s1[1] = p["sigma_solver.S1_DN"];

    const std::string infile = p["solver.INFILE_H5GF"], outfile = p["solver.OUTFILE_H5GF"];
    int n_meshes;
    {
      alps::hdf5::archive in(infile, "r");
      in["/G0/mesh/N"] >> n_meshes;
    }
    if (n_meshes == 3)
      solve<cluster_matsubara_kspace_gf, cluster_itime_kspace_gf>(infile, outfile, NORMAL, s);
    else
      solve<cluster_nambu_matsubara_kspace_gf, cluster_nambu_itime_kspace_gf>(
          infile, outfile, p["sc.NAMBU"].as<bool>() ? NAMBU : DOUBLE_CELL, s);
  } catch (const std::exception& e) {
    std::cerr << "constant_sigma_solver: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
