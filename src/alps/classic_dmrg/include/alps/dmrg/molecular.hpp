// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// WIP/classic-dmrg-chemistry at d4384a7f8, file src/alps/dmrg/molecular.hpp.
// Adaptations: the declaration carries this component's export macro.
// Added here, not in alpsx: the observables and ci_vector requests, the
// spin-resolved two-particle RDM, the observable values and the CI vector
// of the result.
//
// Classic block DMRG for restricted real molecular Hamiltonians: the dmtk
// engine of <alps/dmrg/classic.hpp> on one spinful-fermion site per spatial
// orbital, with every orbital pair coupled, and the one- and two-particle
// reduced density matrices measured on the final state.
#pragma once

#include <alps/dmrg/classic.hpp>
#include <alps/models/terms.hpp>

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace alps::dmrg::classic {

// Restricted, real spatial-orbital integrals, in row-major chemist notation:
// H = sum_pq,a h[p,q] c+[p,a] c[q,a]
//   + 1/2 sum_pqrs,ab eri[p,q,r,s] c+[p,a] c+[r,b] c[s,b] c[q,a]
//   + core_energy. All arrays are owned; orbital indices retain input order.
struct MolecularInput {
    std::size_t orbitals = 0;
    std::size_t alpha_electrons = 0;
    std::size_t beta_electrons = 0;
    std::vector<double> one_body;
    std::vector<double> two_body;
    double core_energy = 0.;
    // Explicit integral pruning; zero preserves every supplied nonzero entry.
    double integral_cutoff = 0.;
    // Further expectation values of the final state, measured with the RDMs:
    // sums of operator products on models::spinful_fermion() sites, site i
    // being orbital i, with real coefficients and at most four sites per
    // product. A product that changes N or Sz has the value zero.
    std::map<std::string, models::TermSum> observables;
    // Also return the final state's CI vector over every determinant of the
    // sector, contracted from the stored block rotations (at most 64
    // orbitals; the work grows with the number of determinants).
    bool ci_vector = false;
};

struct MolecularResult {
    std::size_t orbitals = 0;
    double energy = 0.; // includes core_energy
    // PySCF convention: dm1[p,q] = <c+[q] c[p]>, one array per spin.
    std::vector<double> one_rdm_alpha;
    std::vector<double> one_rdm_beta;
    // dm2[p,q,r,s] = sum_ab <c+[p,a] c+[r,b] c[s,b] c[q,a]>.
    std::vector<double> two_rdm;
    // Its spin blocks, as PySCF's make_rdm12s: aa[p,q,r,s] =
    // <a+[p] a+[r] a[s] a[q]>, ab[p,q,r,s] = <a+[p] b+[r] b[s] a[q]>, bb
    // likewise; two_rdm = aa + ab + ab.transpose(2,3,0,1) + bb to rounding.
    std::vector<double> two_rdm_aa;
    std::vector<double> two_rdm_ab;
    std::vector<double> two_rdm_bb;
    // MolecularInput::observables, by name.
    std::map<std::string, double> observables;
    // With MolecularInput::ci_vector: one entry per determinant of the sector,
    // ordered by alpha then beta occupation word (bit p = orbital p). The
    // amplitudes use the ALPS mode order of pyalps.ed and ALPS::ed: orbitals
    // ascending, each orbital's down mode before its up mode, the state being
    // the product of creators in that order applied to the vacuum. The
    // global sign makes the largest amplitude positive.
    std::vector<std::uint64_t> ci_alpha;
    std::vector<std::uint64_t> ci_beta;
    std::vector<double> ci_coefficients;
    std::vector<Iteration> iterations;
    double final_truncation_error = 0.;
};

// Ground-state classic block DMRG and native 1-/2-RDM measurements, with the
// sweep, truncation, Lanczos, scratch and callback contract of solve(). The
// input fixes N and Sz: options.sector must be empty, num_eigenvalues one and
// arithmetic real. The integrals must have the real-orbital symmetries after
// the explicit cutoff. Fewer than four orbitals are padded with
// one-dimensional neutral sites to meet the engine's four-site layout, without
// adding states. Returns detached RDMs and values, optionally the CI vector,
// but no MPS or checkpoint: nothing can be measured after the call. RDM
// storage and operator preparation scale with the fourth power of the
// orbital count.
ALPS_CLASSIC_DMRG_DECL MolecularResult solve_molecular(MolecularInput const& input, Options const& options = {},
                                                       IterationCallback callback = {});

} // namespace alps::dmrg::classic
