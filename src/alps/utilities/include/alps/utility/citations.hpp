// SPDX-License-Identifier: MIT
#ifndef ALPS_UTILITY_CITATIONS_HPP
#define ALPS_UTILITY_CITATIONS_HPP

#include <alps/utilities_export.h>
#include <iosfwd>
#include <string>

namespace alps {

/// Return the build's compact startup box for a component in CITATIONS.yaml.
/// Includes the preferred framework paper and references for used components.
/// Throws std::invalid_argument for an unknown component; performs no file I/O.
ALPS_UTILITIES_DECL std::string citation_text(const std::string& component = "framework");

/// Print one notice with distinct bibliography entries and all applicable roles.
/// Suppressed only when ALPS_NO_CITATIONS is exactly "1".
ALPS_UTILITIES_DECL void print_citations(std::ostream& out, const std::string& component = "framework");

/// Complete guidance with titles and DOI links, used by --citations.
ALPS_UTILITIES_DECL std::string citation_details(const std::string& component = "framework");
ALPS_UTILITIES_DECL void print_citation_details(std::ostream& out, const std::string& component = "framework");

} // namespace alps
#endif
