// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at c99b526ec^, file
// src/alps/solve/lattice/lattice.cpp (Lattice::to_latticegraph_xml,
// Lattice::to_graph_xml, fnv1a64_hex, alps_graph_for_boundary). Adaptation:
// free functions, vertex <COORDINATE> elements, 17 significant digits
// everywhere, open chain names recognised; see lattice_xml.hpp.
#include <alps/models/lattice_xml.hpp>

#include <cstdint>
#include <cstdio>
#include <sstream>
#include <vector>

namespace alps {
namespace lattice {
namespace {

// Shortest text that reads back to the same double.
std::string number(double value) {
    char buffer[32];
    std::snprintf(buffer, sizeof(buffer), "%.17g", value);
    return buffer;
}

std::string join(std::vector<double> const& values) {
    std::string out;
    for (std::size_t i = 0; i < values.size(); ++i) out += (i ? " " : "") + number(values[i]);
    return out;
}

std::string join(std::vector<int> const& values) {
    std::string out;
    for (std::size_t i = 0; i < values.size(); ++i) out += (i ? " " : "") + std::to_string(values[i]);
    return out;
}

// FNV-1a, 64-bit, lower-case hex: a cross-implementation name, not a
// security primitive.
std::string fnv1a64_hex(std::string const& s) {
    std::uint64_t h = 14695981039346656037ULL;
    for (unsigned char c : s) {
        h ^= c;
        h *= 1099511628211ULL;
    }
    char buffer[32];
    std::snprintf(buffer, sizeof(buffer), "%016llx", static_cast<unsigned long long>(h));
    return buffer;
}

} // namespace

std::string to_latticegraph_xml(Lattice const& lattice, std::string const& name) {
    auto const& cell = lattice.unit_cell();
    int const dim = cell.dimension();
    std::ostringstream os;
    os << "<LATTICES>\n"
       << "<LATTICE name=\"" << name << " base\" dimension=\"" << dim << "\"><BASIS>";
    for (auto const& vector : lattice.basis_vectors()) os << "<VECTOR>" << join(vector) << "</VECTOR>";
    os << "</BASIS></LATTICE>\n"
       << "<UNITCELL name=\"" << name << " cell\" dimension=\"" << dim << "\">\n";
    for (int v = 0; v < cell.num_vertices(); ++v)
        os << "  <VERTEX type=\"" << cell.vertex_types()[static_cast<std::size_t>(v)] << "\"><COORDINATE>"
           << join(cell.vertex_coords()[static_cast<std::size_t>(v)]) << "</COORDINATE></VERTEX>\n";
    for (UnitCell::Edge const& e : cell.edges())
        os << "  <EDGE type=\"" << e.edge_type << "\">"
           << "<SOURCE vertex=\"" << (e.source + 1) << "\" offset=\"" << join(e.source_offset) << "\"/>"
           << "<TARGET vertex=\"" << (e.target + 1) << "\" offset=\"" << join(e.target_offset) << "\"/></EDGE>\n";
    os << "</UNITCELL>\n"
       << "<LATTICEGRAPH name=\"" << name << "\">\n  <FINITELATTICE>\n"
       << "    <LATTICE ref=\"" << name << " base\"/>\n";
    for (int d = 0; d < dim; ++d)
        os << "    <EXTENT dimension=\"" << (d + 1) << "\" size=\"" << lattice.extent()[static_cast<std::size_t>(d)] << "\"/>\n";
    for (int d = 0; d < dim; ++d)
        os << "    <BOUNDARY dimension=\"" << (d + 1) << "\" type=\""
           << to_string(lattice.boundary()[static_cast<std::size_t>(d)]) << "\"/>\n";
    os << "  </FINITELATTICE>\n  <UNITCELL ref=\"" << name << " cell\"/>\n"
       << "</LATTICEGRAPH>\n</LATTICES>\n";
    return os.str();
}

std::string to_graph_xml(Lattice const& lattice) {
    std::ostringstream os;
    os << "<GRAPH vertices=\"" << lattice.num_sites() << "\" edges=\"" << lattice.num_bonds() << "\">";
    for (Site const& s : lattice.sites())
        os << "\n  <VERTEX id=\"" << (s.index + 1) << "\" type=\"" << s.site_type << "\"><COORDINATE>"
           << join(s.coordinate) << "</COORDINATE></VERTEX>";
    for (Bond const& b : lattice.bonds())
        os << "\n  <EDGE source=\"" << (b.source + 1) << "\" target=\"" << (b.target + 1)
           << "\" type=\"" << b.bond_type << "\"/>";
    os << "\n</GRAPH>\n";
    return os.str();
}

std::string custom_lattice_name(Lattice const& lattice) {
    return "alps_custom_lattice_" + fnv1a64_hex(to_latticegraph_xml(lattice, "alps_custom_lattice"));
}

std::string catalog_graph_name(Lattice const& lattice) {
    std::string const& name = lattice.alps_name();
    if (name.empty()) return {};
    auto const& boundary = lattice.boundary();
    auto const all = [&](Boundary b) {
        for (Boundary x : boundary) if (x != b) return false;
        return true;
    };
    if (name == "square lattice" || name == "simple cubic lattice") {
        if (all(Boundary::periodic)) return name;
        if (all(Boundary::open)) return "open " + name;
        return {};
    }
    if (name == "ladder") {
        // Periodic along the legs, always open across the rungs; the
        // "open ladder" sibling opens the legs too.
        if (boundary.size() == 2 && boundary[1] == Boundary::open)
            return boundary[0] == Boundary::periodic ? "ladder" : "open ladder";
        return {};
    }
    if (all(Boundary::periodic)) return name;
    // The chain factories name their open variants themselves ("open chain
    // lattice", "nnn open chain lattice").
    if (all(Boundary::open) && name.find("open ") != std::string::npos) return name;
    return {};
}

} // namespace lattice
} // namespace alps
