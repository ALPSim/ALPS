// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#ifndef ALPS_SOLVERS_HPP
#define ALPS_SOLVERS_HPP

#include <string>

namespace alps {
// <alps/params.hpp> declares the same alias.
namespace params_ns { class params; }
typedef params_ns::params params;

// Link the corresponding ALPS::maxent, ALPS::cthyb, or ALPS::ctint target.
// Solvers run synchronously, write to output_file, and propagate exceptions.
// Each defines the parameters it reads on a copy of `parameters` and throws
// std::invalid_argument, naming them, if any is missing or malformed.
// With an MPI-enabled SDK, the caller must initialize MPI before CT-QMC calls;
// all ranks in MPI_COMM_WORLD must participate.
namespace solvers {
void maxent(params const& parameters, std::string const& output_file);
void cthyb(params const& parameters, std::string const& output_file);
void ctint(params const& parameters, std::string const& output_file);
}
}

#endif
