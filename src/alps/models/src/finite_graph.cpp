// Copyright (C) 2026 ALPS Collaboration. SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/src/lattice/finite_graph.cpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.
#include <alps/models/finite_graph.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace alps {
void FiniteGraph::validate() const {
    auto require = [](bool condition, char const* message) {
        if (!condition) throw std::invalid_argument(message);
    };
    require(!vertices.empty(), "finite graph requires at least one vertex");
    require(vertices.size() <= std::size_t(std::numeric_limits<int>::max()) &&
            edges.size() <= std::size_t(std::numeric_limits<int>::max()),
            "finite graph exceeds the native index range");
    auto embedding_dimension = dimension;
    auto vector = [&](std::vector<double> const& values) {
        if (values.empty()) return;
        if (embedding_dimension == 0) embedding_dimension = values.size();
        require(values.size() == embedding_dimension, "finite graph coordinate dimensions disagree");
        require(std::all_of(values.begin(), values.end(), [](double v) { return std::isfinite(v); }),
                "finite graph coordinates must be finite");
    };
    for (auto const& vertex : vertices) {
        require(vertex.type >= 0, "finite graph vertex types must be nonnegative");
        vector(vertex.coordinate);
    }
    for (auto const& edge : edges) {
        require(edge.source >= 0 && std::size_t(edge.source) < vertices.size() &&
                edge.target >= 0 && std::size_t(edge.target) < vertices.size(),
                "finite graph edge endpoint is outside the vertex array");
        require(edge.type >= 0, "finite graph edge types must be nonnegative");
        vector(edge.displacement);
    }
}
} // namespace alps
