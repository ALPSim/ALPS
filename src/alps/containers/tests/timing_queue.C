/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2003 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/fixed_capacity_deque.h>
#include <boost/timer.hpp>
#include <deque>
#include <iostream>
#include <queue>

const std::size_t n = (2<<24);
const std::size_t m = 16;

int main()
{
  typedef std::queue<int, alps::fixed_capacity_deque<int, m> > Queue0;
  typedef std::queue<int> Queue1;

  boost::timer t0;
  Queue0 queue0;
  int sum0 = 0;
  for (std::size_t i = 0; i < m; ++i) queue0.push(i);
  for (std::size_t i = 0; i < n; ++i) {
    sum0 += queue0.size();
    sum0 += queue0.front();
    sum0 -= queue0.back();
    queue0.pop();
    queue0.push(i + m);
  }
  std::cout << "std::queue with fixed_capacity_deque "
            << t0.elapsed() << " sec\n";

  boost::timer t1;
  Queue1 queue1;
  int sum1 = 0;
  for (std::size_t i = 0; i < m; ++i) queue1.push(i);
  for (std::size_t i = 0; i < n; ++i) {
    sum1 += queue1.size();
    sum1 += queue1.front();
    sum1 -= queue1.back();
    queue1.pop();
    queue1.push(i + m);
  }
  std::cout << "std::queue with std::deque           "
            << t1.elapsed() << " sec\n";

  if (sum0 != sum1) {
    std::cout << "results are inconsistent!\n";
    return -1;
  }

  std::cout << "done.\n";
  return 0;
}
