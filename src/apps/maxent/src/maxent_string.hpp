/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1998-2026 ALPS Collaboration
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Imported from CQMP/Maxent src/maxent_string.hpp at
// e727650a1f0a427e95efebc409d907b4c84011ec; see ../README.md.
#pragma once
#include <algorithm>
#include <cctype>
#include <iomanip>
#include <sstream>
#include <string>

///convert a string to lower case in place (replaces boost::to_lower)
inline void to_lower(std::string &s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
}

///format a double with 17 significant digits, as boost::lexical_cast<std::string> does
inline std::string to_string_exact(double x) {
  std::ostringstream os;
  os << std::setprecision(17) << x;
  return os.str();
}
