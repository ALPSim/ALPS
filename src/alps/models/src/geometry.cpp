// Copyright (C) 2026 by the ALPS collaboration.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/src/lattice/geometry.cpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.

#include <alps/models/geometry.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace alps {

void LatticeGeometry::validate() const
{
    auto require = [](bool condition, char const* message) {
        if (!condition) throw std::invalid_argument(message);
    };
    auto const dim = extent.size();
    // The graph's boundary-crossing record stores two bits per direction.
    require(dim >= 1 && dim <= 4, "lattice geometry requires 1 to 4 dimensions");
    require(boundary.size() == dim && (basis.empty() || basis.size() == dim),
            "lattice geometry extent, boundary and basis dimensions differ");
    require(!vertices.empty(), "lattice geometry requires a nonempty unit cell");
    std::size_t sites = vertices.size();
    std::size_t bonds = edges.size();
    for (std::size_t d = 0; d < dim; ++d) {
        require(extent[d] > 0, "lattice geometry extents must be positive");
        require(boundary[d] == "periodic" || boundary[d] == "open",
                "lattice geometry boundary must be periodic or open");
        require(sites <= std::numeric_limits<int>::max() / std::size_t(extent[d]),
                "lattice geometry exceeds the native site-index range");
        // Bound the expanded graph before allocating it. Open boundaries
        // and self-loops can reduce this upper bound during construction.
        require(bonds <= std::numeric_limits<int>::max() / std::size_t(extent[d]),
                "lattice geometry exceeds the native bond-index range");
        sites *= extent[d];
        bonds *= extent[d];
    }
    auto vector = [&](auto const& v, bool allow_empty = false) {
        require(v.size() == dim || (allow_empty && v.empty()),
                "lattice geometry vector dimension differs from extent");
        for (double value : v)
            require(std::isfinite(value), "lattice geometry coordinates must be finite");
    };
    for (auto const& v : basis) vector(v);
    require(reciprocal_basis.empty() || reciprocal_basis.size() == dim,
            "lattice geometry reciprocal basis dimension differs from extent");
    for (auto const& v : reciprocal_basis) vector(v);
    for (auto const& v : vertices) {
        require(v.type >= 0, "lattice geometry site types must be nonnegative");
        vector(v.coordinate, true); // Catalog vertices may omit the origin.
    }
    for (auto const& edge : edges) {
        require(edge.source >= 0 && std::size_t(edge.source) < vertices.size()
                && edge.target >= 0 && std::size_t(edge.target) < vertices.size(),
                "lattice geometry edge endpoint is outside the unit cell");
        require(edge.type >= 0, "lattice geometry bond types must be nonnegative");
        vector(edge.source_offset);
        vector(edge.target_offset);
        for (std::size_t d = 0; d < dim; ++d) {
            // Cell indices are nonnegative. Check the largest shifted
            // index before the native graph adds offsets and wraps it.
            auto const max_offset = std::numeric_limits<int>::max() - (extent[d] - 1);
            require(edge.source_offset[d] <= max_offset && edge.target_offset[d] <= max_offset,
                    "lattice geometry offset exceeds the native cell-index range");
        }
    }
}

} // namespace alps
