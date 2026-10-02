// SPDX-License-Identifier: MIT
#ifndef ALPS_UTILITY_CITATION_PROVENANCE_HPP
#define ALPS_UTILITY_CITATION_PROVENANCE_HPP

#include <alps/config.h>
#include <string>
#include <vector>
#include <boost/filesystem/path.hpp>

namespace alps {
namespace hdf5 { class archive; }

/// Immutable application-level recommendations from the producing build.
/// An activity labels the contribution; it is not an execution trace.
struct ALPS_DECL citation_snapshot {
  std::string id, component, software_version, catalog_sha256;
  std::string selection_basis, activity, bibliography_cff, notice, request, note, license_note;
  std::vector<std::string> algorithm, implementation, framework;
  bool operator==(const citation_snapshot&) const;
};
typedef std::vector<citation_snapshot> citation_history;

ALPS_DECL citation_snapshot make_citation_snapshot(const std::string& component = "framework",
                                                   const std::string& activity = "calculation");
/// Missing metadata in an old file returns an empty history. Malformed or newer
/// schemas raise an exception, so replacement writers can preserve the old file.
ALPS_DECL citation_history read_citations(hdf5::archive&);
ALPS_DECL citation_history read_citations(const boost::filesystem::path&);
ALPS_DECL void merge_citations(citation_history&, const citation_history&);
/// Append/deduplicate records. Existing complete records are never overwritten.
/// Call inside the same single-writer/file lock boundary as numerical saving.
ALPS_DECL void write_citations(hdf5::archive&, const citation_history&);
/// Replace the citation set when saving a complete, independently owned result.
/// Supply loaded/inherited histories explicitly; filenames do not imply ancestry.
ALPS_DECL void replace_citations(hdf5::archive&, const citation_history&);
ALPS_DECL void write_citations(hdf5::archive&, const std::string& component,
                              const std::string& activity = "calculation");
}
#endif
