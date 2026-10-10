// SPDX-License-Identifier: MIT
#include <alps/utility/citations.hpp>
#include <cstdlib>
#include <ostream>
#include <stdexcept>

namespace {
struct citation_entry {
  const char* component;
  const char* text;
  const char* details;
};

#include <alps/utility/citations_data.inc>
} // namespace

std::string alps::citation_text(const std::string& component) {
  for (const auto& entry : citation_entries)
    if (component == entry.component)
      return entry.text;
  throw std::invalid_argument("Unknown ALPS citation component: " + component);
}

void alps::print_citations(std::ostream& out, const std::string& component) {
  const char* disabled = std::getenv("ALPS_NO_CITATIONS");
  if (disabled && std::string(disabled) == "1") return;
  out << citation_text(component);
}

std::string alps::citation_details(const std::string& component) {
  for (const auto& entry : citation_entries)
    if (component == entry.component)
      return entry.details;
  throw std::invalid_argument("Unknown ALPS citation component: " + component);
}

void alps::print_citation_details(std::ostream& out, const std::string& component) {
  out << citation_details(component);
}
