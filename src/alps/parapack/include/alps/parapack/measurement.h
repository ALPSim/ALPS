/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef PARAPACK_MEASUREMENT_H
#define PARAPACK_MEASUREMENT_H

#include <alps/alea/observables.hpp>
#include <alps/config.h>

namespace alps {

// Add the measurements of a clone to those of the task. With same_weight every
// clone counts as one sample, see merge_random_clone.
ALPS_DECL void merge_clone(alps::alea::observable_set& total, alps::alea::result_set const& clone,
  bool same_weight);

// For averages over disorder: the mean of every observable of the clone is
// one sample of the observable of the same name in total.
ALPS_DECL void merge_random_clone(alps::alea::observable_set& total,
  alps::alea::result_set const& clone);

// For a quantity that does not fluctuate, such as the temperature: it is
// given twice, so that it has an error, of zero.
inline void add_constant(alps::alea::accumulator obs, double val) {
  if (obs.count() < 2) obs << val;
}

} // end namespace alps

#endif // PARAPACK_MEASUREMENT_H
