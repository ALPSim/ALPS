// Copyright (C) 2005-2026 by Emanuel Gull
// ALPS Project: https://alps.comp-phys.org/
// SPDX-License-Identifier: MIT

// Impurity solver for U=0, used by noninteracting_loop_test: without
// interaction, the cluster Green's function is the bare one, G = G0. It is
// called like a real solver, as `free_solver <parameter file> --MU=<mu>`,
// reads G0 from solver.INFILE_H5GF and writes it as G to solver.OUTFILE_H5GF.

#include <iostream>
#include <stdexcept>

#include "alps/hdf5.hpp"
#include "alps/params.hpp"

#include "cluster_transform.h"

int main(int argc, char** argv) {
  try {
    if (argc < 2) throw std::invalid_argument("usage: free_solver <parameter file> [--MU=<mu>]");
    alps::params p(argv[1]);
    p.define<std::string>("solver.INFILE_H5GF", "alps_solver_in.h5gf", "input file");
    p.define<std::string>("solver.OUTFILE_H5GF", "alps_solver_out.h5gf", "output file");

    cluster_matsubara_kspace_gf G_omega;
    cluster_itime_kspace_gf G_tau;
    {
      alps::hdf5::archive in(p["solver.INFILE_H5GF"].as<std::string>(), "r");
      G_omega.load(in, "/G0");
      G_tau.load(in, "/G0_tau");
    }
    alps::hdf5::archive out(p["solver.OUTFILE_H5GF"].as<std::string>(), "w");
    G_omega.save(out, "/G_omega");
    G_tau.save(out, "/G_tau");
  } catch (const std::exception& e) {
    std::cerr << "free_solver: " << e.what() << std::endl;
    return 1;
  }
  return 0;
}
