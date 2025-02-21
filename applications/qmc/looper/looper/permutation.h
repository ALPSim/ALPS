/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1997-2006 by Synge Todo <wistaria@comp-phys.org>
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

#ifndef LOOPER_PERMUTATION_H
#define LOOPER_PERMUTATION_H

#include <algorithm>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <iterator>

namespace looper {

template<class RandomAccessIter, class RandomNumberGenerator>
void random_shuffle(RandomAccessIter first, RandomAccessIter last,
                    RandomNumberGenerator& rng)
{
  using std::iter_swap;
//   boost::variate_generator<RandomNumberGenerator&,
//     boost::uniform_real<> > random(rng, boost::uniform_real<>(0, 1));
  for (typename std::iterator_traits<RandomAccessIter>::difference_type
         n = last - first; n > 1; ++first, --n)
    iter_swap(first, first + (int)(n * rng()));
}


// helper function : guided_sort_01
//
// example
//   input  perm  : 0 1 2 3 4 5 6 7
//          guide : 1 0 0 1 0 1 1 0
//   output perm  : 7 1 2 4 3 5 6 0
//          guide : does not change

template<class RandomAccessIterator0, class RandomAccessIterator1>
void guided_sort_01(RandomAccessIterator0 first_perm,
                    RandomAccessIterator0 last_perm,
                    RandomAccessIterator1 first_guide)
{
  RandomAccessIterator1 last_guide = first_guide + (last_perm - first_perm);
  if (first_perm == last_perm) return;
  --last_perm;
  --last_guide;
  for (;;) {
    if (last_perm - first_perm <= 0) return;
    while (*first_guide == 0) {
      ++first_perm;
      ++first_guide;
      if (first_perm == last_perm) return;
    }
    while (*last_guide == 1) {
      --last_perm;
      --last_guide;
      if (first_perm == last_perm) return;
    }
    std::iter_swap(first_perm, last_perm);
    ++first_perm;
    ++first_guide;
    --last_perm;
    --last_guide;
  }
}


// function : partitioned_random_shuffle
//
// This function generates a permutaion, which satisfies the
// conservation law.  This is used for generating a boundary graph for
// S >= 1 cases.

template<class RandomAccessIter0, class RandomAccessIter1,
         class RandomNumberGenerator>
void partitioned_random_shuffle(RandomAccessIter0 perm_first,
                               RandomAccessIter0 perm_last,
                               RandomAccessIter1 guide0_first,
                               RandomAccessIter1 guide1_first,
                               RandomNumberGenerator& rng)
{
  typedef typename std::iterator_traits<RandomAccessIter0>::difference_type
    diff_type;
  diff_type n = perm_last - perm_first;

  // sort in two sectors (0 and 1) according to values in guide1
  guided_sort_01(perm_first, perm_last, guide1_first);
  diff_type c = std::count(guide1_first, guide1_first + n, 0);

  // shuffle in each sector
  looper::random_shuffle(perm_first, perm_first + c, rng);
  looper::random_shuffle(perm_first + c, perm_last, rng);

  // reorder permutation according to values in guide0
  guided_sort_01(perm_first, perm_last, guide0_first);
}

} // end namespace looper

#endif // LOOPER_PERMUTATION_H
