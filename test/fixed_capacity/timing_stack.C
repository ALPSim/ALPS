/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002-2003 by Synge Todo <wistaria@comp-phys.org>
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
