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

#include <alps/parapack/queue.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/foreach.hpp>
#include <iostream>
#include <queue>

int main()
{
  // check for boost::posix_time

  boost::posix_time::ptime t0 = boost::posix_time::second_clock::local_time();
  std::cerr << "Current time: " << to_simple_string(t0) << std::endl;

  boost::posix_time::ptime t1 = t0 + boost::posix_time::seconds(10);
  std::cerr << "10 seconds after: " << to_simple_string(t1) << std::endl;

  boost::posix_time::ptime t2 = t0 + boost::posix_time::minutes(5);
  std::cerr << "5 minutes after: " << to_simple_string(t2) << std::endl;

  boost::posix_time::ptime t3 = t0 + boost::posix_time::hours(1);
  std::cerr << "One hour after: " << to_simple_string(t3) << std::endl;

  std::cout << "Order comparison: ";
  if (t0 < t1 && t1 < t2 && t2 < t3)
    std::cout << "OK\n";
  else
    std::cout << "error\n";

  std::priority_queue<alps::check_queue_element_t> queue;
  queue.push(alps::next_checkpoint(0, 0, 0, boost::posix_time::seconds(10)));
  queue.push(alps::next_checkpoint(0, 0, 0, boost::posix_time::seconds(2)));
  queue.push(alps::next_checkpoint(0, 0, 0, boost::posix_time::seconds(5)));
  queue.push(alps::next_checkpoint(0, 0, 0, boost::posix_time::seconds(1)));
  std::cerr << to_simple_string(queue.top().time) << std::endl;

  return 0;
}

