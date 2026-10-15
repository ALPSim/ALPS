// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/include/alps/lattice/finite_graph.hpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.
#pragma once
#include <alps/models_export.h>
#include <alps/models/geometry.hpp>
#include <variant>
#include <cstddef>
#include <string>
#include <vector>

namespace alps {

// Explicit finite topology, independent of a unit cell or lattice catalog.
// Vertex and edge order are significant. Parallel edges and their source/target
// orientations are preserved. Empty coordinates/vectors mean unavailable;
// they are not replaced by an invented embedding or inferred displacement.
struct FiniteGraph {
    struct Vertex {
        int type = 0;
        std::vector<double> coordinate;
        bool operator==(Vertex const& o) const { return type == o.type && coordinate == o.coordinate; }
        bool operator!=(Vertex const& o) const { return !(*this == o); }
    };
    struct Edge {
        int source = -1, target = -1;
        int type = 0;
        std::vector<double> displacement;
        bool operator==(Edge const& o) const {
            return source == o.source && target == o.target && type == o.type && displacement == o.displacement;
        }
        bool operator!=(Edge const& o) const { return !(*this == o); }
    };
    std::string name;
    // Zero means no embedding dimension was declared. Any supplied coordinate
    // and displacement vectors must nevertheless agree in dimension.
    std::size_t dimension = 0;
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    void ALPS_MODELS_DECL validate() const;
    bool operator==(FiniteGraph const& o) const {
        return name == o.name && dimension == o.dimension && vertices == o.vertices && edges == o.edges;
    }
    bool operator!=(FiniteGraph const& o) const { return !(*this == o); }
};
// Construction owns exactly one topology representation.
using GraphInput = std::variant<LatticeGeometry, FiniteGraph>;
} // namespace alps
