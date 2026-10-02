// SPDX-License-Identifier: MIT
#include <alps/utility/citation_provenance.hpp>
#include <alps/hdf5.hpp>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <tuple>

namespace {
const std::string root = "/provenance/alps/citations";
struct citation_snapshot_entry {
  const char *calculation_id, *analysis_id, *unspecified_id;
  const char *component, *software_version, *catalog_sha256, *bibliography_cff, *notice;
  const char *request, *note, *license_note, *algorithm, *implementation, *framework;
};
#include <alps/utility/citation_snapshots.inc>

std::vector<std::string> keys(const char* text) {
  std::vector<std::string> result;
  std::istringstream input(text);
  for (std::string key; std::getline(input, key);) result.push_back(key);
  return result;
}
bool digest(const std::string& value) {
  return value.size() == 64 && std::all_of(value.begin(), value.end(), [](char c) {
    return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
  });
}
void validate(const alps::citation_snapshot& s) {
  if (!digest(s.id) || !digest(s.catalog_sha256) || s.component.empty() ||
      s.software_version.empty() || s.bibliography_cff.empty() || s.notice.empty() ||
      s.selection_basis != "application" || s.framework.empty() ||
      (s.activity != "calculation" && s.activity != "analysis" && s.activity != "unspecified"))
    throw std::runtime_error("Invalid ALPS citation snapshot");
}
long long integer_field(alps::hdf5::archive& ar, const std::string& path) {
  if (!ar.is_scalar(path) || !(ar.is_datatype<signed char>(path) || ar.is_datatype<unsigned char>(path) ||
      ar.is_datatype<short>(path) || ar.is_datatype<unsigned short>(path) ||
      ar.is_datatype<int>(path) || ar.is_datatype<unsigned int>(path) ||
      ar.is_datatype<long>(path) || ar.is_datatype<unsigned long>(path) ||
      ar.is_datatype<long long>(path) || ar.is_datatype<unsigned long long>(path)))
    throw std::runtime_error("ALPS citation provenance requires a scalar integer: " + path);
  long long value = 0;
  ar[path] >> value;
  return value;
}
void check_schema(alps::hdf5::archive& ar) {
  if (!ar.is_group(root)) {
    if (ar.is_data(root)) throw std::runtime_error("ALPS citation provenance path is not a group");
    return;
  }
  const auto version = integer_field(ar, root + "/schema_version");
  if (version != 1) throw std::runtime_error("Unsupported ALPS citation provenance schema: " + std::to_string(version));
}
void fields(alps::hdf5::archive& ar, const std::string& path, alps::citation_snapshot& s, bool write) {
#define FIELD(name) if (write) ar.write_utf8(path + "/" #name, s.name); else ar[path + "/" #name] >> s.name
  FIELD(component); FIELD(software_version); FIELD(catalog_sha256);
  FIELD(selection_basis); FIELD(activity); FIELD(bibliography_cff); FIELD(notice);
  FIELD(request); FIELD(note); FIELD(license_note);
#undef FIELD
#define FIELD(name) if (write) ar[path + "/" #name] << s.name; else ar[path + "/" #name] >> s.name
  FIELD(algorithm); FIELD(implementation); FIELD(framework);
#undef FIELD
}
}

bool alps::citation_snapshot::operator==(const citation_snapshot& s) const {
  return std::tie(id, component, software_version, catalog_sha256, selection_basis, activity,
                  bibliography_cff, notice, request, note, license_note, algorithm, implementation, framework)
      == std::tie(s.id, s.component, s.software_version, s.catalog_sha256, s.selection_basis, s.activity,
                  s.bibliography_cff, s.notice, s.request, s.note, s.license_note, s.algorithm, s.implementation, s.framework);
}

alps::citation_snapshot alps::make_citation_snapshot(const std::string& component, const std::string& activity) {
  for (const auto& entry : citation_snapshot_entries) {
    if (component != entry.component) continue;
    citation_snapshot s;
    if (activity == "calculation") s.id = entry.calculation_id;
    else if (activity == "analysis") s.id = entry.analysis_id;
    else if (activity == "unspecified") s.id = entry.unspecified_id;
    else throw std::invalid_argument("Unknown citation activity: " + activity);
    s.component = entry.component; s.software_version = entry.software_version;
    s.catalog_sha256 = entry.catalog_sha256; s.selection_basis = "application";
    s.activity = activity; s.bibliography_cff = entry.bibliography_cff; s.notice = entry.notice;
    s.request = entry.request; s.note = entry.note; s.license_note = entry.license_note;
    s.algorithm = keys(entry.algorithm); s.implementation = keys(entry.implementation);
    s.framework = keys(entry.framework);
    return s;
  }
  throw std::invalid_argument("Unknown ALPS citation component: " + component);
}

alps::citation_history alps::read_citations(hdf5::archive& ar) {
  check_schema(ar);
  citation_history result;
  if (!ar.is_group(root + "/records")) {
    if (ar.is_data(root + "/records")) throw std::runtime_error("ALPS citation records path is not a group");
    return result;
  }
  for (const auto& id : ar.list_children(root + "/records")) {
    const std::string path = root + "/records/" + id;
    if (!ar.is_group(path)) throw std::runtime_error("ALPS citation snapshot path is not a group");
    if (!ar.is_data(path + "/complete")) continue; // interrupted append
    const auto complete = integer_field(ar, path + "/complete");
    if (!complete) continue;
    if (complete != 1) throw std::runtime_error("Invalid ALPS citation completion marker: " + path);
    citation_snapshot s;
    s.id = id;
    fields(ar, path, s, false);
    validate(s);
    result.push_back(s);
  }
  std::sort(result.begin(), result.end(), [](const citation_snapshot& a, const citation_snapshot& b) {
    return a.id < b.id;
  });
  return result;
}
alps::citation_history alps::read_citations(const boost::filesystem::path& filename) {
  hdf5::archive ar(filename.string(), "r");
  return read_citations(ar);
}
void alps::merge_citations(citation_history& target, const citation_history& source) {
  for (const auto& entry : source) {
    validate(entry);
    const auto found = std::find_if(target.begin(), target.end(), [&](const citation_snapshot& s) { return s.id == entry.id; });
    if (found == target.end()) target.push_back(entry);
    else if (!(*found == entry)) throw std::runtime_error("Conflicting ALPS citation snapshot: " + entry.id);
  }
}
void alps::write_citations(hdf5::archive& ar, const citation_history& history) {
  citation_history all = read_citations(ar);
  const auto existing = all;
  merge_citations(all, history); // validate all input/collisions before any mutation
  if (all.empty()) return;
  if (!ar.is_group(root)) ar[root + "/schema_version"] << 1;
  for (auto s : all) {
    if (std::find_if(existing.begin(), existing.end(), [&](const citation_snapshot& old) { return old.id == s.id; }) != existing.end()) continue;
    const std::string path = root + "/records/" + s.id;
    if (ar.is_group(path)) ar.delete_group(path); // only an incomplete record
    fields(ar, path, s, true);
    ar[path + "/complete"] << 1; // commit marker written last
  }
}
void alps::write_citations(hdf5::archive& ar, const std::string& component, const std::string& activity) {
  write_citations(ar, citation_history{make_citation_snapshot(component, activity)});
}
void alps::replace_citations(hdf5::archive& ar, const citation_history& history) {
  auto previous = read_citations(ar); // reject unsupported/damaged old metadata before replacing it
  citation_history validated;
  merge_citations(validated, history);
  merge_citations(previous, validated); // retained identities must keep the same payload
  if (validated.empty()) {
    if (ar.is_group(root)) ar.delete_group(root);
    return;
  }
  // Retain complete records in place; a checkpoint should not rewrite immutable
  // payloads merely because the caller owns the full result's citation set.
  if (ar.is_group(root + "/records")) {
    for (const auto& id : ar.list_children(root + "/records")) {
      if (std::none_of(validated.begin(), validated.end(), [&](const citation_snapshot& s) { return s.id == id; }))
        ar.delete_group(root + "/records/" + id);
    }
  }
  write_citations(ar, validated);
}
