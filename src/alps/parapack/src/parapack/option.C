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
#include <alps/parapack/option.h>
#include <boost/lexical_cast.hpp>
#include <iostream>

namespace alps {
namespace parapack {

namespace pt = boost::posix_time;

option::option(int argc, char** argv, bool eval)
  : desc("Allowed options"), for_evaluate(eval), show_help(false),
    show_license(false), time_limit(pt::pos_infin),
    check_interval(pt::millisec(100)), checkpoint_interval(pt::seconds(3600)),
    report_interval(pt::seconds(600)), vmusage_interval(pt::pos_infin),
    use_termfile(false), auto_evaluate(true), evaluate_only(false),
    dump_format(dump_format::hdf5), dump_policy(dump_policy::RunningOnly),
    task_range(), write_xml(false),
    use_mpi(false), default_total_threads(true), auto_total_threads(false),
    num_total_threads(1), threads_per_clone(1), jobfiles(), valid(true) {
  desc.flag("help,h", "produce help message")
    .flag("license,l", "print license conditions")
    .value<std::string>("dump-format",
     "format for dumping info, parameter, and measurements [hdf5 (default), xdr]")
    .value<std::string>("task-range",
     "specify range of task indices to be processed, e.g. [2:5]")
    .flag("write-xml", "write results to XML files")
    .list("input-file", "input master XML files")
    .positional("input-file", -1);
  if (!for_evaluate) {
    desc.flag("auto-evaluate", "evaluate observables upon halting [default = true]")
      .flag("check-parameter", "perform parameter checking")
      .value<int>("check-interval",
       "time between internal status check [unit = millisec; default = 100ms]")
      .value<int>("checkpoint-interval",
       "time between checkpointing [unit = sec; default = 3600s]")
      .value<std::string>("dump-policy",
       "policy for dumping user checkpoint data [running (default), never, all]")
      .flag("enable-termination-file", "enable termination file support (*.term)")
      .flag("evaluate", "evaluation mode")
      .flag("mpi", "run in parallel using MPI")
      .value<int>("Nmin", "obsolete")
      .value<int>("Nmax", "obsolete")
      .flag("no-evaluate", "prevent evaluating observables upon halting")
      .value<int>("report-interval",
       "time between progress report of clones [unit = sec; default = 600s]")
      .value<int>("vmusage-interval",
       "time between virtual memory usage report of processes [unit = sec; default = none]")
      .value<int>("time-limit,T",
       "time limit for the simulation [unit = sec; default = no time limit]")
      .value<int>("Tmin", "obsolete")
      .value<int>("Tmax", "obsolete")
      .value<int>("threads-per-clone,p",
       "number of threads for each clone [default = 1]")
      .value<std::string>("total-threads,r",
       "total number of threads [integer or 'auto'; default = total number of processes]");
  }

  try {
    desc.parse(argc, argv);
  }
  catch (std::exception& e) {
    valid = false;
    std::cerr << e.what() << std::endl;
    return;
  }

  if (desc.count("help"))
    show_help = true;
  if (desc.count("license"))
    show_license = true;
  if (desc.count("dump-format")) {
    std::string value = desc.get<std::string>("dump-format");
    if (value == "hdf5")
      dump_format = dump_format::hdf5;
    else if (value == "xdr")
      dump_format = dump_format::xdr;
    else {
      valid = false;
      return;
    }
  }
  if (desc.count("task-range"))
    task_range = task_range_t(desc.get<std::string>("task-range"));
  if (desc.count("write-xml"))
    write_xml = true;
  if (desc.count("input-file"))
    jobfiles = desc.values("input-file");
  if (!for_evaluate) {
    if (desc.count("auto-evaluate"))
      auto_evaluate = true;
    if (desc.count("no-evaluate"))
      auto_evaluate = false;
    if (desc.count("check-interval"))
      check_interval = pt::millisec(desc.get<int>("check-interval"));
    if (desc.count("checkpoint-interval"))
      checkpoint_interval = pt::seconds(desc.get<int>("checkpoint-interval"));
    if (desc.count("enable-termination-file"))
      use_termfile = true;
    if (desc.count("dump-policy")) {
      std::string value = desc.get<std::string>("dump-policy");
      if (value == "never")
        dump_policy = dump_policy::Never;
      else if (value == "running")
        dump_policy = dump_policy::RunningOnly;
      else if (value == "all")
        dump_policy = dump_policy::All;
      else {
        valid = false;
        return;
      }
    }
    if (desc.count("report-interval"))
      report_interval = pt::seconds(desc.get<int>("report-interval"));
    if (desc.count("vmusage-interval"))
      vmusage_interval = pt::seconds(desc.get<int>("vmusage-interval"));
    if (desc.count("mpi"))
      use_mpi = true;
    if (desc.count("evaluate"))
      evaluate_only = true;
    if (desc.count("time-limit"))
      time_limit = pt::seconds(desc.get<int>("time-limit"));
    if (desc.count("threads-per-clone"))
      threads_per_clone = desc.get<int>("threads-per-clone");
    if (desc.count("total-threads")) {
      default_total_threads = false;
      std::string const value = desc.get<std::string>("total-threads");
      if (value == "auto") {
        auto_total_threads = true;
      } else {
        try {
          num_total_threads = boost::lexical_cast<int>(value);
        }
        catch (boost::bad_lexical_cast const&) {
          valid = false;
          std::cerr << "the argument ('" << value << "') for option '--total-threads' is invalid"
                    << std::endl;
          return;
        }
      }
    }
  }
}

void option::print(std::ostream& os) const { desc.print(os); }

void option::print_summary(std::ostream& os, std::string const& prefix) const {
  if (!for_evaluate) {
    os << prefix << "auto evaluation = " << (auto_evaluate ? "yes" : "no") << std::endl;
    os << prefix << "time limit = ";
    if (!time_limit.is_special())
      os << time_limit.total_seconds() << " seconds\n";
    else
      os << "unlimited\n";
    os << prefix << "interval between checkpointing  = "
       << checkpoint_interval.total_seconds() << " seconds\n";
    os << prefix << "interval between progress report = "
       << report_interval.total_seconds() << " seconds\n";
    os << prefix << "interval between vmusage report = ";
    if (!vmusage_interval.is_special())
      os << vmusage_interval.total_seconds() << " seconds\n";
    else
      os << "infinity\n";
  }
  os << prefix << "task range = ";
  if (task_range.valid())
    os << task_range << std::endl;
  else
    os << "all\n";
  os << prefix << "worker dump format = " << dump_format::to_string(dump_format) << std::endl;
  if (!for_evaluate) {
    os << prefix << "worker dump policy = " << dump_policy::to_string(dump_policy) << std::endl;
  }
}

} // end namespace parapack
} // end namespace alps
