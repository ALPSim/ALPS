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

#ifndef PARAPACK_CLONE_TIMER_H
#define PARAPACK_CLONE_TIMER_H

#include <alps/config.h>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace alps {

class clone_timer {
public:
  typedef uint64_t loops_t;
  typedef boost::posix_time::ptime time_t;
  typedef boost::posix_time::time_duration duration_t;

  clone_timer(duration_t const& check_interval, double progress = 0) {
    interval_ = check_interval;
    reset(progress);
  }

  void reset(double progress = 0) {
    start_time_ = current_time();
    start_progress_ = progress;
    next_check_ = start_time_ + interval_;
  }

  static time_t current_time() { return boost::posix_time::microsec_clock::local_time(); }

  loops_t next_loops(loops_t loops) { return next_loops(loops, current_time()); }
  loops_t next_loops(loops_t loops, time_t const& current) {
    if (current > next_check_) {
      loops /= 2;
      if (loops == 0) loops = 1;
    } else if (current + interval_/2 < next_check_) {
      loops *= 2;
    }
    next_check_ = current + interval_;
    return loops;
  }

private:
  duration_t interval_;
  time_t start_time_;
  double start_progress_;
  time_t next_check_;
};

} // end namespace alps

#endif // PARAPACK_CLONE_TIMER_H
