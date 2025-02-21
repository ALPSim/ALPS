/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2010 by Synge Todo <wistaria@comp-phys.org>
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
