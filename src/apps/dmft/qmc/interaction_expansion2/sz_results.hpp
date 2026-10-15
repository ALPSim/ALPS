// Copyright (C) 2026 ALPS Collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT

#ifndef DMFT_QMC_INTERACTION_EXPANSION2_SZ_RESULTS_HPP
#define DMFT_QMC_INTERACTION_EXPANSION2_SZ_RESULTS_HPP

#include <alps/accumulators.hpp>

/// The Sz_i, Sz2_i and Sz0_Szi accumulators of sites 0 to n_site-1 measure
/// the value times the sign. ALPS 3.0 measured them as signed observables,
/// whose results are divided by the mean sign. Divide each of these results
/// that has measurements by the "Sign" result, so they keep that meaning.
void divide_sz_results_by_sign(alps::accumulators::result_set & results, unsigned int n_site);

#endif
