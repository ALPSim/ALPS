/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/
#include <alps/parapack/measurement.h>

#include <string>

namespace alps {

void merge_clone(alps::alea::observable_set& total, alps::alea::result_set const& clone,
                 bool same_weight) {
  if (same_weight) {
    merge_random_clone(total, clone);
  } else {
    total.merge(clone);
  }
}

void merge_random_clone(alps::alea::observable_set& total, alps::alea::result_set const& clone) {
  for (std::string const& name : clone.names()) {
    alps::alea::observable_result const& obs = clone[name];
    if (obs.is_histogram() || obs.count() == 0) continue;
    if (!total.has(name)) total << alps::alea::observable(name);
    if (obs.is_scalar())
      total[name] << obs.mean();
    else
      total[name] << obs.mean_vector();
  }
}

} // end namespace alps
