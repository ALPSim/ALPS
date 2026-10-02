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

#ifndef PARAPACK_PERMUTATION_H
#define PARAPACK_PERMUTATION_H

#include <algorithm>
#include <boost/random/uniform_real.hpp>
#include <iterator>

namespace alps {

template<class RandomAccessIter, class RandomNumberGenerator>
void random_shuffle(RandomAccessIter first, RandomAccessIter last, RandomNumberGenerator& rng) {
  using std::iter_swap;
  for (typename std::iterator_traits<RandomAccessIter>::difference_type
         n = last - first; n > 1; ++first, --n)
    iter_swap(first, first + (int)(n * rng()));
}

} // end namespace alps

#endif // PARAPACK_PERMUTATION_H
