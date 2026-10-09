// SPDX-License-Identifier: MIT
#ifndef ALPS_UTILITY_CLI_HPP
#define ALPS_UTILITY_CLI_HPP

#include <alps/config.h>
#include <functional>
#include <string>

namespace alps {

/// Own MPI initialization only when the caller has not already initialized it.
/// Used by informational CLI queries and legacy drivers that print before their
/// solver initializes MPI. Does not select a simulation's execution mode.
class ALPS_DECL cli_mpi_guard {
public:
  cli_mpi_guard(int& argc, char**& argv);
  ~cli_mpi_guard();
  cli_mpi_guard(const cli_mpi_guard&) = delete;
  cli_mpi_guard& operator=(const cli_mpi_guard&) = delete;
private:
  bool owns_mpi_;
};

/// True for serial execution or MPI rank zero, independent of scheduler mode.
ALPS_DECL bool cli_is_master();

/// Handle a standalone --help/-h, --license/-l, or --citations query on rank zero.
/// An optional --mpi is accepted. Multiple queries or calculation arguments are
/// rejected. Returns false for a calculation invocation (including after --).
ALPS_DECL bool handle_cli_information(int argc, char** argv,
    const std::string& component, const std::function<void()>& print_help);

} // namespace alps
#endif
