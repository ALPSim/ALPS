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

#include <alps/fixed_capacity_vector.h>
#include <boost/timer.hpp>
#include <deque>
#include <iostream>
#include <stack>

const std::size_t n = (2<<24);
const std::size_t m = 16;

int main()
{
  typedef std::stack<int, alps::fixed_capacity_vector<int, m> > Stack0;
  typedef std::stack<int> Stack1;

  boost::timer t0;
  Stack0 stack0;
  int sum0 = 0;
  for (std::size_t i = 0; i < m; ++i) stack0.push(i);
  for (std::size_t i = 0; i < n; ++i) {
    sum0 -= stack0.size();
    sum0 += stack0.top();
    stack0.pop();
    stack0.push(i + m);
  }
  std::cout << "std::stack with fixed_capacity_vector "
            << t0.elapsed() << " sec\n";

  boost::timer t1;
  Stack1 stack1;
  int sum1 = 0;
  for (std::size_t i = 0; i < m; ++i) stack1.push(i);
  for (std::size_t i = 0; i < n; ++i) {
    sum1 -= stack1.size();
    sum1 += stack1.top();
    stack1.pop();
    stack1.push(i + m);
  }
  std::cout << "std::stack with std::deque            "
            << t1.elapsed() << " sec\n";

  if (sum0 != sum1) {
    std::cout << "results are inconsistent!\n";
    return -1;
  }

  std::cout << "done.\n";
  return 0;
}
