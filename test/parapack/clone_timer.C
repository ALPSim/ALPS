/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2009 by Synge Todo <wistaria@comp-phys.org>
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

#include <alps/parapack/clone_timer.h>
#include <iostream>

int main()
{
#ifndef BOOST_NO_EXCEPTIONS
try {
#endif

  double progress = 0;
  alps::clone_timer timer(boost::posix_time::milliseconds(100), progress);

  alps::clone_timer::loops_t loops = 1024;
  alps::clone_timer::time_t start = timer.current_time();
  alps::clone_timer::time_t end = start + boost::posix_time::seconds(1);
  std::cerr << "start time = " << start << std::endl;

  while (true) {
    std::cerr << "loops = " << loops << std::endl;
    for (alps::clone_timer::loops_t i = 0; i < loops; ++i) {
      timer.current_time();
    }
    alps::clone_timer::time_t current = timer.current_time();
    std::cerr << "current time = " << current << " (" << (current - start) << ")\n";
    if (current > end) break;
    loops = timer.next_loops(loops, current);
  }
  end += boost::posix_time::seconds(1);
  while (true) {
    std::cerr << "loops = " << loops << std::endl;
    for (alps::clone_timer::loops_t i = 0; i < loops; ++i) {
      timer.current_time();
      timer.current_time();
      timer.current_time();
      timer.current_time();
    }
    alps::clone_timer::time_t current = timer.current_time();
    std::cerr << "current time = " << current << " (" << (current - start) << ")\n";
    if (current > end) break;
    loops = timer.next_loops(loops, current);
  }

#ifndef BOOST_NO_EXCEPTIONS
}
catch (std::exception& exp) {
  std::cerr << exp.what() << std::endl;
  std::abort();
}
#endif
  return 0;
}
