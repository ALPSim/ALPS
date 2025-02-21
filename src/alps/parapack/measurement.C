/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2008 by Synge Todo <wistaria@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#include "measurement.h"

namespace alps {

namespace {

template<typename U, typename T>
bool merge_obs(alps::ObservableSet& total, alps::Observable const* obs) {
  if (dynamic_cast<T const*>(obs)) {
    if (dynamic_cast<T const*>(obs)->count()) {
      std::string name = obs->name();
      if (!total.has(name)) {
        total.addObservable(U(name));
        total[name].reset(true);
      }
      total[name] << dynamic_cast<T const*>(obs)->mean();
    }
    return true;
  }
  return false;
}

} // end namespace

void merge_clone(alps::ObservableSet& total, alps::ObservableSet const& clone, bool same_weight) {
  if (same_weight) {
    merge_random_clone(total, clone);
  } else {
    total << clone;
  }
}

void merge_random_clone(ObservableSet& total, ObservableSet const& clone) {
  for (std::map<std::string, Observable*>::const_iterator itr = clone.begin();
       itr != clone.end(); ++itr) {
    if (itr->second) {
      merge_obs<SimpleRealObservable, SimpleRealObservable>(total, itr->second) ||
      merge_obs<RealObservable, RealObservable>(total, itr->second) ||
      merge_obs<RealObservable, RealTimeSeriesObservable>(total, itr->second) ||
      merge_obs<RealObservable, SignedObservable<SimpleRealObservable> >(total, itr->second) ||
      merge_obs<RealObservable, SignedObservable<RealObservable> >(total, itr->second) ||
      merge_obs<RealObservable, SignedObservable<RealTimeSeriesObservable> >(total, itr->second) ||
      merge_obs<RealObservable, RealObsevaluator>(total, itr->second) ||
      merge_obs<RealObservable, AbstractSignedObservable<RealObsevaluator> >(total, itr->second) ||

      merge_obs<RealVectorObservable, SimpleRealVectorObservable>(total, itr->second) ||
      merge_obs<RealVectorObservable, RealVectorObservable>(total, itr->second) ||
      merge_obs<RealVectorObservable, RealVectorTimeSeriesObservable>(total, itr->second) ||
      merge_obs<RealVectorObservable, SignedObservable<SimpleRealVectorObservable>
        >(total, itr->second) ||
      merge_obs<RealVectorObservable, SignedObservable<RealVectorObservable>
        >(total, itr->second) ||
      merge_obs<RealVectorObservable, SignedObservable<RealVectorTimeSeriesObservable>
        >(total, itr->second) ||
      merge_obs<RealVectorObservable, RealVectorObsevaluator>(total, itr->second) ||
      merge_obs<RealVectorObservable, AbstractSignedObservable<RealVectorObsevaluator>
        >(total, itr->second);
    }
  }
}

} // end namespace alps
