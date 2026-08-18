/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2013 by Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "logger.h"
#include <boost/lexical_cast.hpp>
#include <boost/date_time.hpp>

namespace alps {

std::string logger::header() {
  return std::string("[") + to_simple_string(boost::posix_time::second_clock::local_time()) + "]: ";
}
std::string logger::task(alps::tid_t tid) {
  return std::string("task[") + boost::lexical_cast<std::string>(tid+1) + ']';
}
std::string logger::clone(alps::tid_t tid, alps::cid_t cid) {
  return std::string("clone[") + boost::lexical_cast<std::string>(tid+1) + ',' +
    boost::lexical_cast<std::string>(cid+1) + ']';
}
std::string logger::group(alps::process_group g) {
  return std::string("processgroup[") + boost::lexical_cast<std::string>(g.group_id+1) + ']';
}
std::string logger::group(alps::thread_group g) {
  return std::string("threadgroup[") + boost::lexical_cast<std::string>(g.group_id+1) + ']';
}

std::string logger::usage(alps::vmusage_type const& u) {
  return std::string("Process ID = ") +
    boost::lexical_cast<std::string>(u.find("Pid")->second) + ", " +
    "VmPeak = " + boost::lexical_cast<std::string>(u.find("VmPeak")->second) + " [kB], " +
    "VmSize = " + boost::lexical_cast<std::string>(u.find("VmSize")->second) + " [kB], " +
    "VmHWM = " + boost::lexical_cast<std::string>(u.find("VmHWM")->second) + " [kB], " +
    "VmRSS = " + boost::lexical_cast<std::string>(u.find("VmRSS")->second) + " [kB]";
}

} // namespace alps
