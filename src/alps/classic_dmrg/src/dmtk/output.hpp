// Copyright (C) 2026 by the ALPS collaboration.
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/solvers/tensor_network/src/classic_dmrg/dmtk/output.hpp, the Eigen-based
// rework of src/apps/dmrg/dmrg/dmtk in this tree (Adrian Feiguin's dmtk).
// Adaptation: "dmtk/" include paths; dmtk::half_integer and dmtk/support.hpp
// replace the ALPS model and utility headers.

#pragma once

#include <iostream>

namespace dmtk {

// Only classic DMRG's calling thread is redirected. Unowned CLI execution
// retains its original streams; no process-wide stream buffer is replaced.
inline thread_local std::ostream* output_stream = nullptr;
inline thread_local std::ostream* error_stream = nullptr;

inline std::ostream& out() noexcept {
  return output_stream ? *output_stream : std::cout;
}

inline std::ostream& err() noexcept {
  return error_stream ? *error_stream : std::cerr;
}

} // namespace dmtk
