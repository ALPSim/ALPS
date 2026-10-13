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

#include <alps/alea.h>

namespace alps {

ALPS_DECL void merge_clone(alps::ObservableSet& total, alps::ObservableSet const& clone,
  bool same_weight);

ALPS_DECL void merge_random_clone(alps::ObservableSet& total, alps::ObservableSet const& clone);

template<typename T>
void add_constant(Observable& obs, T const& val) {
  if (dynamic_cast<SimpleRealObservable*>(&obs) &&
      dynamic_cast<SimpleRealObservable*>(&obs)->count() < 2) obs << val;
}

} // end namespace alps

#endif // PARAPACK_MEASUREMENT_H
