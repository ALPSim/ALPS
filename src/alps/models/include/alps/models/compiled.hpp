// Copyright (C) 2026 ALPS Collaboration
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/include/alps/model/compiled.hpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports.
#pragma once

#include <complex>
#include <map>
#include <string>
#include <vector>

namespace alps::models {

// One located matrix, in row-major order (row=output, column=input).
// Its dimension is supplied by CompiledOperator::site_dims[site].
struct LocalFactor {
    int site = 0;
    std::vector<std::complex<double>> matrix;
    bool operator==(LocalFactor const& o) const { return site == o.site && matrix == o.matrix; }
    bool operator!=(LocalFactor const& o) const { return !(*this == o); }
};

// A coefficient times an ordered product of local factors and parity strings.
// An empty factor list is a scalar identity. The complete record keeps factor
// locations and matrices together; there are no parallel arrays to synchronize.
struct LocalTerm {
    std::complex<double> coefficient{1.};
    std::vector<LocalFactor> factors;
    std::vector<int> jordan_wigner;
    bool operator==(LocalTerm const& o) const {
        return coefficient == o.coefficient && factors == o.factors && jordan_wigner == o.jordan_wigner;
    }
    bool operator!=(LocalTerm const& o) const { return !(*this == o); }
};

// Owned numeric operator data, independent of lattice catalogs, measurements,
// solver selection and execution. It can describe a non-Hermitian observable.
// Consumers validate their own scientific restrictions without modifying it.
struct CompiledOperator {
    std::vector<int> site_dims;
    std::vector<LocalTerm> terms;
    bool operator==(CompiledOperator const& o) const { return site_dims == o.site_dims && terms == o.terms; }
    bool operator!=(CompiledOperator const& o) const { return !(*this == o); }
};

// Connectivity in an already tiled finite lattice. Parallel edges are retained.
struct ModelEdge {
    int src;
    int dst;
    int bond_type;
    bool operator==(ModelEdge const& o) const { return src == o.src && dst == o.dst && bond_type == o.bond_type; }
    bool operator!=(ModelEdge const& o) const { return !(*this == o); }
};

// Model metadata composed with the operator. ED can consume op directly;
// lattice samplers additionally require the geometry and local-space metadata.
struct CompiledModel {
    CompiledOperator op;
    std::vector<int> site_types;
    std::vector<ModelEdge> edges;
    // Physical basis eigenvalues in local-state order, indexed by site type.
    // These retain their meaning when measurement operators are overridden.
    std::map<int, std::vector<double>> spin_basis;
    std::map<int, std::vector<double>> number_basis;
    // Requested diagonal measurement operators, indexed by site type/name.
    std::map<int, std::map<std::string, std::vector<double>>> diag_observables;
};

} // namespace alps::models
