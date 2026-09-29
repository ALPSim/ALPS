// SPDX-License-Identifier: MIT
#ifndef ALPS_UTILITY_CITATIONS_HPP
#define ALPS_UTILITY_CITATIONS_HPP

#include <alps/config.h>
#include <iosfwd>
#include <string>

namespace alps {

/// Return the build's citation notice for a component in CITATIONS.yaml.
/// Includes the preferred framework paper and references for used components.
/// Throws std::invalid_argument for an unknown component; performs no file I/O.
ALPS_DECL std::string citation_text(const std::string& component = "framework");

/// Print one notice with distinct bibliography entries and all applicable roles.
ALPS_DECL void print_citations(std::ostream& out, const std::string& component = "framework");

} // namespace alps
#endif
