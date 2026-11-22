// Copyright (C) 2005-2026 by Emanuel Gull
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

// Impurity solver for U=0, used by noninteracting_loop_test: without
// interaction, the cluster Green's function is the bare one, G = G0. It is
// called like a real solver, as `free_solver <parameter file> --MU=<mu>`,
// reads G0 from solver.INFILE_H5GF and writes it as G to solver.OUTFILE_H5GF.

#include <iostream>
#include <string>
#include <stdexcept>

#include "alps/hdf5.hpp"
#include "alps/params.hpp"

#include "cluster_transform.h"

namespace {
// Copies G0 to G for one framework's Green's function types.
template<class GFW, class GFT> void copy_g0(const std::string& infile, const std::string& outfile) {
  GFW G_omega;
  GFT G_tau;
  {
    alps::hdf5::archive in(infile, "r");
    G_omega.load(in, "/G0");
    G_tau.load(in, "/G0_tau");
  }
  alps::hdf5::archive out(outfile, "w");
  G_omega.save(out, "/G_omega");
  G_tau.save(out, "/G_tau");
}
}

int main(int argc, char** argv) {
  try {
    if (argc < 2) throw std::invalid_argument("usage: free_solver <parameter file> [--MU=<mu>]");
    // the parameter file and the --NAME=value overrides dca passes on
    alps::params p(argc, (const char**)argv);
    p.define<std::string>("solver.INFILE_H5GF", "alps_solver_in.h5gf", "input file");
    p.define<std::string>("solver.OUTFILE_H5GF", "alps_solver_out.h5gf", "output file");
    const std::string infile = p["solver.INFILE_H5GF"], outfile = p["solver.OUTFILE_H5GF"];

    // The normal-state G0 has the meshes (frequency, K, flavor); the double
    // unit cell and Nambu G0 have a fourth, a 2x2 index.
    int n_meshes;
    {
      alps::hdf5::archive in(infile, "r");
      in["/G0/mesh/N"] >> n_meshes;
    }
    if (n_meshes == 3)
      copy_g0<cluster_matsubara_kspace_gf, cluster_itime_kspace_gf>(infile, outfile);
    else
      copy_g0<cluster_nambu_matsubara_kspace_gf, cluster_nambu_itime_kspace_gf>(infile, outfile);
  } catch (const std::exception& e) {
    std::cerr << "free_solver: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
