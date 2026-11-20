/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2003 by Matthias Troyer <troyer@itp.phys.ethz.ch>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#include <alps/scheduler/options.h>
#include <alps/utilities/copyright.hpp>
#include <alps/cli.hpp>
#include <alps/cli/options.hpp>
#include <boost/limits.hpp>
#include <boost/throw_exception.hpp>
#include <iostream>
#include <stdexcept>

namespace alps {
namespace scheduler {

namespace {

// The options every scheduler program takes; Options adds the job file.
cli::options declare_options(NoJobfileOptions& o)
{
  cli::options desc("Allowed options");
  desc.flag("help", "produce help message")
    .flag("license,l", "print license conditions")
    .flag("mpi", "run in parallel using MPI")
    .value("checkpoint-time", o.checkpoint_time, 1800, "time between checkpoints")
    .value("Tmin", o.min_check_time, 60, "minimum time between checks whether a simulation is finished")
    .value("Tmax", o.max_check_time, 900, "maximum time between checks whether a simulation is finished")
    .value("time-limit,T", o.time_limit, 0, "time limit for the simulation")
    .value("Nmin", o.min_cpus, 1, "minimum number of CPUs per simulation")
    .value("Nmax", o.max_cpus, (std::numeric_limits<int>::max)(), "maximum number of CPUs per simulation")
    .flag("write-xml", "write results to XML files");
  return desc;
}

}

NoJobfileOptions::NoJobfileOptions()
  : min_check_time(60),
    max_check_time(900),
    checkpoint_time(1800),
    min_cpus(1),
    max_cpus(1),
    time_limit(0.),
    use_mpi(false),
    valid(true), // shall we really run?
    write_xml(false)
{
}

NoJobfileOptions::NoJobfileOptions(int argc, char** argv)
  : NoJobfileOptions(argc, argv, "scheduler") {}

NoJobfileOptions::NoJobfileOptions(int argc, char** argv, const std::string& citation_component)
  : programname(std::string(argv[0])),
    use_mpi(false),
    valid(true), // shall we really run?
    write_xml(false)
{
  if (argc) {
  
  cli::options desc = declare_options(*this);
  if (alps::handle_cli_information(argc, argv, citation_component, [&] { std::cout << desc << "\n"; })) {
    valid = false;
    return;
  }
  desc.parse(argc, argv);

  if (desc.count("mpi")) {
    use_mpi = true;
  }

#ifdef ALPS_HAVE_HDF5
  if (desc.count("write-xml"))
    write_xml = true;
#else
  write_xml=true;
#endif

  if(min_cpus>max_cpus)
    boost::throw_exception(std::runtime_error("Minimum number of CPUs larger than maximum number of CPU"));
  if(min_check_time>max_check_time)
    boost::throw_exception(std::runtime_error("Minimum time between checks larger than maximum time"));
  }
}

Options::Options()
  : NoJobfileOptions()
{
}

Options::Options(int argc, char** argv)
  : Options(argc, argv, "scheduler") {}

Options::Options(int argc, char** argv, const std::string& citation_component)
{
  programname = std::string(argv[0]);
  valid = true;
  if (argc) {
  std::string filename;
  
  cli::options desc = declare_options(*this);
  desc.value("input-file", filename, "input file").positional("input-file");
  if (alps::handle_cli_information(argc, argv, citation_component, [&] { std::cout << desc << "\n"; })) {
    valid = false;
    return;
  }
  desc.parse(argc, argv);

  if (desc.count("mpi")) {
    use_mpi = true;
  }

#ifdef ALPS_HAVE_HDF5
  if (desc.count("write-xml"))
    write_xml = true;
#else
  write_xml=true;
#endif

  if (!filename.empty())
    jobfilename=std::filesystem::path(filename);
  else
    boost::throw_exception(std::runtime_error("No job file specified"));
    
  if(min_cpus>max_cpus)
    boost::throw_exception(std::runtime_error("Minimum number of CPUs larger than maximum number of CPU"));
  if(min_check_time>max_check_time)
    boost::throw_exception(std::runtime_error("Minimum time between checks larger than maximum time"));
  }
}
    
} // namespace scheduler
} // namespace alps
