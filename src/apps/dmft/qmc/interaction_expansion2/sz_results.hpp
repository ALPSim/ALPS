// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT

#ifndef DMFT_QMC_INTERACTION_EXPANSION2_SZ_RESULTS_HPP
#define DMFT_QMC_INTERACTION_EXPANSION2_SZ_RESULTS_HPP

#include <alps/accumulators.hpp>

/// Every observable of the Green's function, density and Sz measurements is
/// accumulated as value times the configuration's sign. Divide each of these
/// results that has measurements by the "Sign" result (with jackknife), so
/// they are the physical expectation values <O> = <O sign> / <sign>.
void divide_signed_results_by_sign(alps::accumulators::result_set & results,
                                   unsigned int n_site, unsigned int n_flavors);

#endif
