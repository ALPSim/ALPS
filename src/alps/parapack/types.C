/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1997-2011 by Synge Todo <wistaria@comp-phys.org>
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

#include "types.h"

namespace alps {

std::string dump_format::to_string(dump_format_t format) {
  switch (format) {
  case hdf5:
    return "hdf5";
  case xdr:
    return "xdr";
  default:
    ;
  }
  return "undefined";
}

std::string dump_policy::to_string(dump_policy_t policy) {
  switch (policy) {
  case All:
    return "all";
  case RunningOnly:
    return "running workers only";
  case Never:
    return "never";
  default:
    ;
  }
  return "undefined";
}

// task_status_t task_status::status(double progress, double max_work, bool on_memory) {
//   if (on_memory)
//     if (progress < 1)
//       return Running;
//     else
//       return (progress < max_work) ? Continuing : Idling;
//   else
//     if (progress < 1)
//       return Suspended;
//     else
//       return (progress < max_work) ? Finished : Completed;
// }

task_status_t task_status::status(std::string const& str) {
  if (str == "new")
    return NotStarted;
  else if (str == "running")
    return Suspended;
  else if (str == "finished")
    return Finished;
  else if (str == "completed")
    return Completed;
  else
    boost::throw_exception(std::runtime_error("invalid status string"));
  return Undefined;
}

std::string task_status::to_string(task_status_t status) {
  switch (status) {
  case Ready:
    return "new";
  case Running:
    return "running";
  case Continuing:
    return "finished";
  case Idling:
    return "completed";
  case NotStarted:
    return "new";
  case Suspended:
    return "running";
  case Finished:
    return "finished";
  case Completed:
    return "completed";
  default:
    boost::throw_exception(std::runtime_error("invalid status"));
  }
  return "undefined";
}

clone_status_t clone_status::status(double progress) {
  return (progress < 1) ? Running : Finished;
}

clone_status_t clone_status::status(double progress, bool on_memory) {
  if (on_memory)
    return (progress < 1) ? Running : Idling;
  else
    return (progress < 1) ? Suspended : Finished;
}

clone_status_t clone_status::status(std::string const& str) {
  if (str == "running")
    return clone_status::Suspended;
  else if (str == "finished")
    return clone_status::Finished;
  else
    boost::throw_exception(std::runtime_error("invalid status string"));
  return Undefined;
}

std::string clone_status::to_string(clone_status_t status) {
  switch (status) {
  case Running:
    return "running";
  case Finished:
    return "finished";
  default:
    boost::throw_exception(std::runtime_error("invalid status"));
  }
  return "undefined";
}

} // end namespace alps
