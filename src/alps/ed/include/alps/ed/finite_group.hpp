// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/solvers/diag/include/alps/diag/finite_group.hpp.
// Adaptation: C++17 (no std::span), <alps/ed/...> include names, ALPS_ED_DECL
// exports, namespace alps::ed.
#pragma once

#include <alps/ed_export.h>

#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace alps::ed {

// Immutable finite site-permutation group with one unitary irreducible
// representation. Permutations map source to destination; row-major matrices
// obey D(g h)=D(g)D(h), with (g h)[i]=g[h[i]]. Closure checks every generator
// relation, unitarity and character norm <chi,chi>=1 (absolute tolerance 1e-10).
// Identity is first. The site action must be faithful; local rotations,
// projective and antiunitary representations are outside this contract.
class ALPS_ED_DECL FiniteGroupIrrep {
public:
    FiniteGroupIrrep(std::vector<std::vector<std::size_t>> generators,
                     std::vector<std::vector<std::complex<double>>> generator_matrices,
                     std::size_t dimension, std::string label, std::size_t max_order = 4096);
    std::size_t order() const;
    std::size_t num_sites() const;
    std::size_t dimension() const;
    std::string const& label() const;
    std::vector<std::vector<std::size_t>> const& generators() const;
    /// Every group element as a site permutation, identity first.
    std::vector<std::vector<std::size_t>> const& permutations() const;
    /// The row-major dimension x dimension matrix of every element.
    std::vector<std::vector<std::complex<double>>> const& matrices() const;
    std::vector<std::complex<double>> const& characters() const;

private:
    struct Data;
    std::shared_ptr<const Data> data_;
};

// D_n on 3<=n<=2048 ring sites: rotation i->i+1, reflection i->-i.
// k in [0,n/2]. k=0 and even-n k=n/2 are 1D (reflection=+/-1);
// interior k gives a real 2D irrep and requires reflection=+1.
ALPS_ED_DECL FiniteGroupIrrep dihedral_irrep(std::size_t n, std::size_t angular_index = 0, int reflection = 1);
// Require a complete, nonduplicate family of the same permutation group.
// Element order may differ; labels must be unique and nonempty.
ALPS_ED_DECL void validate_irrep_family(std::vector<FiniteGroupIrrep> const&);

} // namespace alps::ed
