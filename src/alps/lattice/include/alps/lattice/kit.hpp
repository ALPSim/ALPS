// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// The bridge from ALPS::models lattices (alps::lattice::Lattice) to the
// legacy lattice library: a kit lattice becomes either the catalog
// LATTICEGRAPH it is identical to, with the L/W/H and bc parameters that
// select it, or a custom <LATTICEGRAPH> document carrying the kit's exact
// structure. Either form drives graph_helper in process or, written to a
// file, the XML applications through LATTICE_LIBRARY. Header-only; the
// documents come from <alps/models/lattice_xml.hpp>, so users link
// ALPS::models besides ALPS::alps.
//
// Follows the to_alps_params design of skilledwolf/alpsx (branch
// agent/v3-execution-and-api-audit at c99b526ec^, src/alps/solve/lattice/
// lattice.cpp), with two changes: a named catalog graph is used only when the
// catalog has it and expands to the same sites, types, coordinates and
// ordered bonds as the kit lattice (alpsx trusted the name), and the custom
// document is handed back as text instead of an in-process registry.

#ifndef ALPS_LATTICE_KIT_HPP
#define ALPS_LATTICE_KIT_HPP

#include <alps/lattice/graph_helper.h>
#include <alps/lattice/latticelibrary.h>
#include <alps/models/lattice.hpp>
#include <alps/models/lattice_xml.hpp>
#include <alps/parameter.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace alps {
namespace lattice {

struct LatticeResolution {
  /// LATTICE and, on the catalog route, L/W/H, bc1..bcD and the catalog's
  /// own LATTICE_LIBRARY when one was given.
  Parameters parameters;
  /// The LATTICEGRAPH name selected.
  std::string graph;
  /// The <LATTICES> document on the custom route; empty on the catalog route.
  std::string library_xml;
  bool from_catalog() const { return library_xml.empty(); }
};

/// Why a legacy graph differs from a kit lattice, "" when it does not: site
/// and bond counts, site types, coordinates within `tolerance`, and the
/// multiset of ordered (source, target, type) bonds.
template <class G>
std::string structure_difference(graph_helper<G> const& graph, Lattice const& lattice, double tolerance = 1e-12) {
  if (static_cast<int>(graph.num_sites()) != lattice.num_sites())
    return "site count " + std::to_string(graph.num_sites()) + " != " + std::to_string(lattice.num_sites());
  if (static_cast<int>(graph.num_bonds()) != lattice.num_bonds())
    return "bond count " + std::to_string(graph.num_bonds()) + " != " + std::to_string(lattice.num_bonds());
  for (auto it = graph.sites().first; it != graph.sites().second; ++it) {
    auto const index = static_cast<std::size_t>(*it);
    Site const& site = lattice.sites()[index];
    if (static_cast<int>(graph.site_type(*it)) != site.site_type) return "site type at site " + std::to_string(index);
    auto const& coordinate = graph.coordinate(*it);
    if (coordinate.size() != site.coordinate.size()) return "coordinate dimension at site " + std::to_string(index);
    for (std::size_t d = 0; d < coordinate.size(); ++d)
      if (!(std::abs(coordinate[d] - site.coordinate[d]) <= tolerance))
        return "coordinate at site " + std::to_string(index);
  }
  using Edge = std::tuple<int, int, int>;
  std::vector<Edge> legacy, kit;
  for (auto it = graph.bonds().first; it != graph.bonds().second; ++it)
    legacy.emplace_back(static_cast<int>(graph.source(*it)), static_cast<int>(graph.target(*it)),
                        static_cast<int>(graph.bond_type(*it)));
  for (Bond const& bond : lattice.bonds()) kit.emplace_back(bond.source, bond.target, bond.bond_type);
  std::sort(legacy.begin(), legacy.end());
  std::sort(kit.begin(), kit.end());
  if (legacy != kit) return "bonds (source, target, type) differ";
  return {};
}

namespace kit_detail {
inline Parameters catalog_parameters(Lattice const& lattice, std::string const& graph, Parameters const& catalog) {
  Parameters p;
  if (catalog.defined("LATTICE_LIBRARY")) p["LATTICE_LIBRARY"] = catalog["LATTICE_LIBRARY"];
  p["LATTICE"] = graph;
  static char const* const extents[] = {"L", "W", "H"};
  for (int d = 0; d < lattice.dimension() && d < 3; ++d)
    p[extents[d]] = lattice.extent()[static_cast<std::size_t>(d)];
  for (int d = 0; d < lattice.dimension(); ++d)
    p["bc" + std::to_string(d + 1)] = to_string(lattice.boundary()[static_cast<std::size_t>(d)]);
  return p;
}
} // namespace kit_detail

/// Always the custom document: the kit's structure verbatim, for in-process
/// use where no catalog name is wanted.
inline LatticeResolution custom_lattice(Lattice const& lattice) {
  LatticeResolution r;
  r.graph = custom_lattice_name(lattice);
  r.library_xml = to_latticegraph_xml(lattice, r.graph);
  r.parameters["LATTICE"] = r.graph;
  return r;
}

/// The catalog route when the lattice library selected by `catalog`
/// (LATTICE_LIBRARY, else lattices.xml) has catalog_graph_name(lattice) and
/// that graph, expanded with the lattice's extent, has the kit's structure;
/// the custom document otherwise: for a lattice without catalog name, a name
/// the catalog lacks (the kit's "kagome lattice" against the catalog's
/// "Kagome lattice", the rectangular and non-cubic 3-D graphs), and a graph
/// that differs (the untyped honeycomb and kagome vertices, boundaries the
/// catalog does not bake).
inline LatticeResolution resolve_lattice(Lattice const& lattice, Parameters const& catalog = Parameters()) {
  std::string const name = catalog_graph_name(lattice);
  if (!name.empty() && lattice.dimension() <= 3) {
    LatticeLibrary const library(catalog);
    // graph_helper accepts a LATTICEGRAPH or a flat GRAPH of that name.
    if (library.has_lattice(name) || library.has_graph(name)) {
      Parameters const p = kit_detail::catalog_parameters(lattice, name, catalog);
      graph_helper<> const expanded(p);
      if (structure_difference(expanded, lattice).empty()) {
        LatticeResolution r;
        r.graph = name;
        r.parameters = p;
        return r;
      }
    }
  }
  return custom_lattice(lattice);
}

/// The legacy graph of a resolution, built in process: the custom document is
/// parsed from memory, the catalog route reads the lattice library. `extra`
/// parameters are merged first; the resolution's keys take precedence.
template <class G = typename graph_helper<>::graph_type>
std::unique_ptr<graph_helper<G>> make_graph_helper(LatticeResolution const& resolution,
                                                   Parameters const& extra = Parameters()) {
  Parameters p = extra;
  p << resolution.parameters;
  if (resolution.from_catalog()) return std::make_unique<graph_helper<G>>(p);
  std::istringstream in(resolution.library_xml);
  return std::make_unique<graph_helper<G>>(in, p);
}

/// Write the custom document to `file` (replacing it) and return the
/// resolution's parameters with LATTICE_LIBRARY set to that file's absolute
/// path, as an XML application expects them. The catalog route writes
/// nothing and returns the parameters unchanged.
inline Parameters write_lattice_library(LatticeResolution const& resolution, std::filesystem::path const& file) {
  if (resolution.from_catalog()) return resolution.parameters;
  std::ofstream out(file, std::ios::binary | std::ios::trunc);
  if (!out) throw std::runtime_error("cannot write the lattice library " + file.string());
  out << resolution.library_xml;
  out.close();
  if (!out) throw std::runtime_error("cannot write the lattice library " + file.string());
  Parameters p = resolution.parameters;
  p["LATTICE_LIBRARY"] = std::filesystem::absolute(file).string();
  return p;
}

} // namespace lattice
} // namespace alps

#endif // ALPS_LATTICE_KIT_HPP
