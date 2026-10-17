// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/include/alps/model/builtin.hpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.
//
// Built-in model factories. Same models, term forms (Sz-eigenoperator
// bond expressions for QMC compatibility), coupling-name metadata,
// conserved QNs and defaults as builtin.py. Coefficients are numeric
// (bound) — the symbolic layer stays Python-side.

#pragma once

#include <alps/models_export.h>
#include <alps/models/hamiltonian.hpp>

namespace alps {
namespace models {

// -- spin models --
ALPS_MODELS_DECL Hamiltonian heisenberg(lattice::Lattice const& lat, double J = 1.0, double h = 0.0,
                       double S = 0.5, std::optional<int> bond_type = std::nullopt);
ALPS_MODELS_DECL Hamiltonian ising(lattice::Lattice const& lat, double J = 1.0, double h = 0.0,
                  double S = 0.5, std::optional<int> bond_type = std::nullopt);
ALPS_MODELS_DECL Hamiltonian xxz(lattice::Lattice const& lat, double Jxy = 1.0, double Jz = 1.0,
                double h = 0.0, double S = 0.5,
                std::optional<int> bond_type = std::nullopt);
ALPS_MODELS_DECL Hamiltonian xy(lattice::Lattice const& lat, double J = 1.0, double S = 0.5,
               std::optional<int> bond_type = std::nullopt);
ALPS_MODELS_DECL Hamiltonian transverse_field_ising(lattice::Lattice const& lat, double J = 1.0,
                                   double Gamma = 1.0, double h = 0.0, double S = 0.5,
                                   std::optional<int> bond_type = std::nullopt);
ALPS_MODELS_DECL Hamiltonian heisenberg_d(lattice::Lattice const& lat, double J = 1.0, double D = 0.0,
                         double h = 0.0, double S = 1.0,
                         std::optional<int> bond_type = std::nullopt);
ALPS_MODELS_DECL Hamiltonian biquadratic(lattice::Lattice const& lat, double J = 1.0, double K = 0.0,
                        double S = 1.0, std::optional<int> bond_type = std::nullopt);
// Bond type 0 = nearest (J1), 1 = next-nearest (J2); use with nnn_chain /
// frustrated_square.
ALPS_MODELS_DECL Hamiltonian j1j2(lattice::Lattice const& lat, double J1 = 1.0, double J2 = 0.5,
                 double S = 0.5);

// -- fermionic models --
ALPS_MODELS_DECL Hamiltonian hubbard(lattice::Lattice const& lat, double t = 1.0, double U = 1.0,
                    double mu = 0.0, std::optional<int> bond_type = std::nullopt);
ALPS_MODELS_DECL Hamiltonian spinless_fermion_model(lattice::Lattice const& lat, double t = 1.0,
                                   double mu = 0.0, double V = 0.0,
                                   std::optional<int> bond_type = std::nullopt);
ALPS_MODELS_DECL Hamiltonian t_j(lattice::Lattice const& lat, double t = 1.0, double J = 0.4,
                std::optional<int> bond_type = std::nullopt);
// Two site types: 0 = itinerant SpinfulFermion, 1 = localized spin-1/2.
// Hopping on bond type 0, Kondo exchange on bond type 1.
ALPS_MODELS_DECL Hamiltonian kondo_lattice(lattice::Lattice const& lat, double t = 1.0, double J = 1.0,
                          double mu = 0.0);

// -- bosonic models --
ALPS_MODELS_DECL Hamiltonian hardcore_boson_model(lattice::Lattice const& lat, double t = 1.0,
                                 double mu = 0.0, double V = 0.0,
                                 std::optional<int> bond_type = std::nullopt);
ALPS_MODELS_DECL Hamiltonian bose_hubbard(lattice::Lattice const& lat, double t = 1.0, double U = 1.0,
                         double mu = 0.0, double V = 0.0, int n_max = 3,
                         std::optional<int> bond_type = std::nullopt);

// -- molecular Hamiltonians --
// H = sum_{pq,s} h_pq c+_{ps} c_{qs} + 1/2 sum_{pqrs,st} (pq|rs) c+_{ps} c+_{rt} c_{st} c_{qs}
// on one spinful_fermion() site per spatial orbital of an open chain (strict
// convention, N and Sz conserved). h1 is the n x n one-electron matrix and eri
// the n^4 two-electron integrals in chemist's notation, both row-major.
// Integrals that are zero or smaller in magnitude than `tol` are dropped. The
// terms, their order and their coefficients are those of
// qdk_chemistry_alps.abinitio.build_abinitio_hamiltonian, as one added sum.
ALPS_MODELS_DECL Hamiltonian ab_initio(std::vector<double> const& h1, std::vector<double> const& eri,
                                       int n_orbitals, double tol = 1e-12);

} // namespace models
} // namespace alps
