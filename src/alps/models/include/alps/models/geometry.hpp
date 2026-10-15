// Copyright (C) 2026 by the ALPS collaboration.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/include/alps/lattice/geometry.hpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.
#pragma once

#include <alps/models_export.h>
#include <string>
#include <vector>

namespace alps {

// A finite repeated cell, in numeric coordinates. This owns its inputs and
// requires no catalog, parameter expressions, files, or process registry.
// Basis vectors and fractional cell coordinates retain their input precision.
struct LatticeGeometry {
    struct Vertex {
        int type = 0;
        std::vector<double> coordinate;
    };
    struct Edge {
        int source;
        int target;
        int type;
        std::vector<int> source_offset;
        std::vector<int> target_offset;
    };

    std::vector<int> extent;
    std::vector<std::string> boundary;
    // Empty preserves an unspecified physical embedding (e.g. unit-cell input).
    // It must not be replaced by an invented Cartesian basis.
    std::vector<std::vector<double>> basis;
    // Empty means momentum components are expressed in the lattice basis.
    // When supplied, vectors include the factor 2*pi, as in the catalog.
    std::vector<std::vector<double>> reciprocal_basis;
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;

    void ALPS_MODELS_DECL validate() const;

};

} // namespace alps
