// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#ifndef ALPS_SOLVERS_HPP
#define ALPS_SOLVERS_HPP

#include <string>

namespace alps {
class params;

// Link the corresponding ALPS::maxent, ALPS::cthyb, or ALPS::ctint target.
// Solvers run synchronously, write to output_file, and propagate exceptions.
// With an MPI-enabled SDK, the caller must initialize MPI before CT-QMC calls;
// all ranks in MPI_COMM_WORLD must participate.
namespace solvers {
void maxent(params const& parameters, std::string const& output_file);
void cthyb(params const& parameters, std::string const& output_file);
void ctint(params const& parameters, std::string const& output_file);
}
}

#endif
