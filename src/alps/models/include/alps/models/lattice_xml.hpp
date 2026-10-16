// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at c99b526ec^ (the last revision before
// the fork switched its lattice handoff to TOML), the to_latticegraph_xml,
// to_graph_xml, fnv1a64_hex and alps_graph_for_boundary parts of
// src/alps/solve/lattice/lattice.cpp. Adaptation: free functions on
// alps::lattice::Lattice; unit-cell vertices carry their fractional
// <COORDINATE> (alpsx emitted none, so multi-vertex cells lost their
// coordinates); every value is written with 17 significant digits; the
// catalog graph name also recognises the open chain variants the factories
// pre-resolve; the in-process registry and alps::Parameters bridge live in
// <alps/lattice/kit.hpp> of the legacy lattice module.
//
// XML descriptions of kit lattices in the dialect of the ALPS lattice
// library (lattices.xml), so that programmatic lattices can drive the XML
// applications and the legacy graph_helper.

#pragma once

#include <alps/models_export.h>
#include <alps/models/lattice.hpp>

#include <string>

namespace alps {
namespace lattice {

// The complete <LATTICES> document of one lattice: a <LATTICE name="<name>
// base"> with the basis vectors, a <UNITCELL name="<name> cell"> whose typed
// vertices carry their fractional <COORDINATE> and whose edges carry their
// cell offsets, and a <LATTICEGRAPH name="<name>"> whose <FINITELATTICE>
// bakes the extent and the per-dimension boundary. The legacy expansion
// enumerates sites and bonds in the kit's order (cells with the last
// dimension fastest, one bond per cell and unit-cell edge) and reproduces the
// kit's coordinates; the one divergence is a periodic extent of 1, where the
// legacy expansion keeps the self-loop the kit drops.
ALPS_MODELS_DECL std::string to_latticegraph_xml(Lattice const& lattice,
                                                 std::string const& name = "alps_custom_lattice");

// The flat <GRAPH vertices=".." edges=".."> of the tiled cluster: every site
// with its 1-based id, type and real-space <COORDINATE>, every bond with its
// 1-based source and target and its type. No unit-cell metadata.
ALPS_MODELS_DECL std::string to_graph_xml(Lattice const& lattice);

// "alps_custom_lattice_<16 hex digits>": FNV-1a over the document emitted
// under the placeholder name, so equal geometries share one name and a
// parameter sweep over geometries can keep several in one library.
ALPS_MODELS_DECL std::string custom_lattice_name(Lattice const& lattice);

// The catalog LATTICEGRAPH that bakes this lattice's boundary: alps_name()
// when every direction is periodic or the factory already named the open
// variant (the chains), the "open ..." sibling of the square, simple cubic
// and ladder graphs, "" when the catalog has no graph for this boundary or
// the lattice is custom. Neither existence in a catalog nor structural
// identity is checked here; <alps/lattice/kit.hpp> does both.
ALPS_MODELS_DECL std::string catalog_graph_name(Lattice const& lattice);

} // namespace lattice
} // namespace alps
