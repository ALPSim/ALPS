// Copyright (C) 2026 ALPS collaboration. SPDX-License-Identifier: MIT
#include <alps/params.hpp>
#include <alps/solvers.hpp>

#include <stdexcept>
#include <string>

// Calling the public entry point pulls the complete static solver into this
// executable. Link only ALPS::maxent to verify its exported dependency closure.
int main() {
  alps::params parameters;
  parameters.supply("OMEGA_MAX", "4").supply("BETA", "2").supply("NDAT", "0")
            .supply("NFREQ", "20").supply("NORM", "1").supply("N_ALPHA", "4")
            .supply("ALPHA_MIN", "0.1").supply("ALPHA_MAX", "10");
  try {
    alps::solvers::maxent(parameters, "unused.out.h5");
  } catch (std::invalid_argument const& error) {
    if (std::string(error.what()) != "NDAT too small") return 1;
  }
  // A missing parameter is an error naming it.
  alps::params incomplete;
  incomplete.supply("BETA", "2");
  try {
    alps::solvers::maxent(incomplete, "unused.out.h5");
  } catch (std::invalid_argument const& error) {
    return std::string(error.what()).find("NDAT") != std::string::npos ? 0 : 3;
  }
  return 2;
}
