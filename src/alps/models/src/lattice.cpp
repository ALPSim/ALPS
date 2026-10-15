// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file src/model/src/construction/lattice.cpp.
// Adaptation: C++17 (no std::span, defaulted comparisons or std::erase_if), <alps/models/...> include names, ALPS_MODELS_DECL exports;
// legacy Parameters bridge, TOML emitter and lattice registry removed.
//
// Implementation of alps::lattice (see lattice.hpp). Semantics mirror the
// Python alps.lattice module line-for-line where behavior is observable;
// the golden-fixture parity suite (tests/lattice_kit_parity.cpp)
// pins the equivalence.

#include <alps/models/lattice.hpp>

#include <cmath>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <utility>


namespace alps {
namespace lattice {

namespace {

// Python's modulo: non-negative result for a positive modulus.
inline int py_mod(int x, int n) { return ((x % n) + n) % n; }

std::vector<Boundary> uniform_bc(Boundary bc, int dim) {
    return std::vector<Boundary>(static_cast<std::size_t>(dim), bc);
}

} // namespace

// -- Boundary -----------------------------------------------------------------

std::string to_string(Boundary b) {
    return b == Boundary::periodic ? "periodic" : "open";
}

Boundary boundary_from_string(std::string const& s) {
    if (s == "periodic") return Boundary::periodic;
    if (s == "open") return Boundary::open;
    throw std::invalid_argument("alps::lattice: unknown boundary '" + s +
                                "' (expected 'periodic' or 'open')");
}

// -- UnitCell -----------------------------------------------------------------

UnitCell::UnitCell(int dimension) : dimension_(dimension) {
    if (dimension < 1)
        throw std::invalid_argument("alps::lattice::UnitCell: dimension must be >= 1");
}

int UnitCell::add_vertex(std::vector<double> frac_coord, int vertex_type) {
    if (static_cast<int>(frac_coord.size()) != dimension_)
        throw std::invalid_argument(
            "alps::lattice::UnitCell::add_vertex: coordinate length != dimension");
    vertex_coords_.push_back(std::move(frac_coord));
    vertex_types_.push_back(vertex_type);
    return num_vertices() - 1;
}

void UnitCell::add_edge(int source, int target, int edge_type,
                        std::vector<int> source_offset,
                        std::vector<int> target_offset) {
    auto const check_vertex = [this](int v, char const* what) {
        if (v < 0 || v >= num_vertices())
            throw std::out_of_range(std::string("alps::lattice::UnitCell::add_edge: ") +
                                    what + " vertex out of range");
    };
    check_vertex(source, "source");
    check_vertex(target, "target");
    auto const norm_offset = [this](std::vector<int>& o, char const* what) {
        if (o.empty()) o.assign(static_cast<std::size_t>(dimension_), 0);
        if (static_cast<int>(o.size()) != dimension_)
            throw std::invalid_argument(std::string("alps::lattice::UnitCell::add_edge: ") +
                                        what + " offset length != dimension");
    };
    norm_offset(source_offset, "source");
    norm_offset(target_offset, "target");
    edges_.push_back(Edge{source, target, edge_type,
                          std::move(source_offset), std::move(target_offset)});
}

Lattice UnitCell::tile(std::vector<int> const& extent,
                       std::vector<std::vector<double>> basis_vectors,
                       std::vector<Boundary> bc) const {
    if (static_cast<int>(extent.size()) != dimension_)
        throw std::invalid_argument(
            "alps::lattice::UnitCell::tile: extent length != dimension");
    if (basis_vectors.empty()) {
        basis_vectors.assign(static_cast<std::size_t>(dimension_),
                             std::vector<double>(static_cast<std::size_t>(dimension_), 0.0));
        for (int d = 0; d < dimension_; ++d) basis_vectors[d][d] = 1.0;
    }
    if (bc.empty()) bc = uniform_bc(Boundary::periodic, dimension_);
    if (static_cast<int>(bc.size()) != dimension_)
        throw std::invalid_argument(
            "alps::lattice::UnitCell::tile: bc length != dimension");
    return Lattice::build(*this, std::move(basis_vectors), extent, std::move(bc));
}

Lattice UnitCell::tile(std::vector<int> const& extent,
                       std::vector<std::vector<double>> basis_vectors,
                       Boundary bc_uniform) const {
    return tile(extent, std::move(basis_vectors), uniform_bc(bc_uniform, dimension_));
}

// -- Lattice ------------------------------------------------------------------

Lattice::Lattice(std::vector<Site> sites, std::vector<Bond> bonds,
                 std::vector<std::vector<double>> basis_vectors,
                 std::vector<int> extent, std::vector<Boundary> boundary,
                 UnitCell unit_cell, std::string alps_name)
    : sites_(std::move(sites)), bonds_(std::move(bonds)),
      basis_vectors_(std::move(basis_vectors)), extent_(std::move(extent)),
      boundary_(std::move(boundary)), unit_cell_(std::move(unit_cell)),
      alps_name_(std::move(alps_name)) {
    neighbors_.resize(sites_.size());
    for (Bond const& b : bonds_) {
        neighbors_[static_cast<std::size_t>(b.source)].push_back(b.target);
        neighbors_[static_cast<std::size_t>(b.target)].push_back(b.source);
    }
}

Lattice Lattice::build(UnitCell unit_cell,
                       std::vector<std::vector<double>> basis_vectors,
                       std::vector<int> extent,
                       std::vector<Boundary> boundary,
                       std::string alps_name) {
    int const dim = unit_cell.dimension();
    int const nv = unit_cell.num_vertices();
    if (static_cast<int>(extent.size()) != dim || static_cast<int>(boundary.size()) != dim
        || static_cast<int>(basis_vectors.size()) != dim)
        throw std::invalid_argument("alps::lattice: extent, boundary and basis must match the cell dimension");
    std::size_t const coordinate_dim = basis_vectors.front().size();
    if (coordinate_dim == 0)
        throw std::invalid_argument("alps::lattice: basis vectors must have at least one component");
    for (auto const& vector : basis_vectors)
        if (vector.size() != coordinate_dim)
            throw std::invalid_argument("alps::lattice: basis vectors must have consistent widths");
    if (nv == 0)
        throw std::invalid_argument("alps::lattice: unit cell has no vertices");
    for (int n : extent)
        if (n < 1) throw std::invalid_argument("alps::lattice: extent entries must be >= 1");

    // Cell enumeration: lexicographic with the LAST dimension fastest —
    // exactly itertools.product(range(n_0), ..., range(n_{d-1})).
    std::vector<std::vector<int>> cell_indices;
    {
        std::vector<int> idx(static_cast<std::size_t>(dim), 0);
        while (true) {
            cell_indices.push_back(idx);
            int d = dim - 1;
            for (; d >= 0; --d) {
                if (++idx[static_cast<std::size_t>(d)] < extent[static_cast<std::size_t>(d)]) break;
                idx[static_cast<std::size_t>(d)] = 0;
            }
            if (d < 0) break;
        }
    }
    // cell index tuple -> flat id (== its position in cell_indices).
    std::map<std::vector<int>, int> cell_id;
    for (std::size_t flat = 0; flat < cell_indices.size(); ++flat)
        cell_id[cell_indices[flat]] = static_cast<int>(flat);

    // -- sites --
    std::vector<Site> sites;
    sites.reserve(cell_indices.size() * static_cast<std::size_t>(nv));
    for (std::size_t ci = 0; ci < cell_indices.size(); ++ci) {
        std::vector<int> const& cell = cell_indices[ci];
        for (int vi = 0; vi < nv; ++vi) {
            std::vector<double> coord(coordinate_dim, 0.0);
            for (int d = 0; d < dim; ++d) {
                double const w = static_cast<double>(cell[static_cast<std::size_t>(d)]) +
                                 unit_cell.vertex_coords()[static_cast<std::size_t>(vi)]
                                                          [static_cast<std::size_t>(d)];
                for (std::size_t c = 0; c < coordinate_dim; ++c)
                    coord[c] += w * basis_vectors[static_cast<std::size_t>(d)][c];
            }
            Site s;
            s.index = static_cast<int>(ci) * nv + vi;
            s.site_type = unit_cell.vertex_types()[static_cast<std::size_t>(vi)];
            s.sublattice = vi;
            s.coordinate = std::move(coord);
            sites.push_back(std::move(s));
        }
    }

    // -- bonds --
    // One bond per (cell, unit-cell edge) instance — NO unordered-pair
    // deduplication, matching the ALPS FINITELATTICE expansion
    // bond-for-bond (convention decision 2026-07-02: a periodic extent-2
    // direction yields BOTH wrap images; any offset o with 2o ≡ 0 (mod L)
    // doubles its pairs). Self-bonds are still dropped — L=1-periodic
    // stays outside the legacy bridge's validated regime.
    std::vector<Bond> bonds;
    for (std::vector<int> const& cell : cell_indices) {
        for (UnitCell::Edge const& edge : unit_cell.edges()) {
            std::vector<int> s_cell(static_cast<std::size_t>(dim));
            std::vector<int> t_cell(static_cast<std::size_t>(dim));
            for (int d = 0; d < dim; ++d) {
                s_cell[static_cast<std::size_t>(d)] =
                    cell[static_cast<std::size_t>(d)] + edge.source_offset[static_cast<std::size_t>(d)];
                t_cell[static_cast<std::size_t>(d)] =
                    cell[static_cast<std::size_t>(d)] + edge.target_offset[static_cast<std::size_t>(d)];
            }
            std::vector<int> s_wrapped = s_cell;
            std::vector<int> t_wrapped = t_cell;
            bool skip = false;
            for (int d = 0; d < dim && !skip; ++d) {
                int const n = extent[static_cast<std::size_t>(d)];
                if (boundary[static_cast<std::size_t>(d)] == Boundary::periodic) {
                    s_wrapped[static_cast<std::size_t>(d)] = py_mod(s_cell[static_cast<std::size_t>(d)], n);
                    t_wrapped[static_cast<std::size_t>(d)] = py_mod(t_cell[static_cast<std::size_t>(d)], n);
                } else {
                    if (s_cell[static_cast<std::size_t>(d)] < 0 ||
                        s_cell[static_cast<std::size_t>(d)] >= n) { skip = true; break; }
                    if (t_cell[static_cast<std::size_t>(d)] < 0 ||
                        t_cell[static_cast<std::size_t>(d)] >= n) { skip = true; break; }
                }
            }
            if (skip) continue;

            auto const s_it = cell_id.find(s_wrapped);
            auto const t_it = cell_id.find(t_wrapped);
            if (s_it == cell_id.end() || t_it == cell_id.end()) continue;

            int const s_idx = s_it->second * nv + edge.source;
            int const t_idx = t_it->second * nv + edge.target;
            if (s_idx == t_idx) continue;

            Bond b;
            b.source = s_idx;
            b.target = t_idx;
            b.bond_type = edge.edge_type;
            b.vector.resize(coordinate_dim);
            for (std::size_t d = 0; d < coordinate_dim; ++d)
                b.vector[d] = sites[static_cast<std::size_t>(t_idx)].coordinate[d] -
                              sites[static_cast<std::size_t>(s_idx)].coordinate[d];
            b.wraps = (s_wrapped != s_cell) || (t_wrapped != t_cell);
            bonds.push_back(std::move(b));
        }
    }

    return Lattice(std::move(sites), std::move(bonds), std::move(basis_vectors),
                   std::move(extent), std::move(boundary), std::move(unit_cell),
                   std::move(alps_name));
}

std::vector<Site> Lattice::sites_by_type(int site_type) const {
    std::vector<Site> out;
    for (Site const& s : sites_) if (s.site_type == site_type) out.push_back(s);
    return out;
}

std::vector<Bond> Lattice::bonds_by_type(int bond_type) const {
    std::vector<Bond> out;
    for (Bond const& b : bonds_) if (b.bond_type == bond_type) out.push_back(b);
    return out;
}

std::vector<int> Lattice::periodic_directions() const {
    std::vector<int> out;
    for (int d = 0; d < dimension(); ++d)
        if (boundary_[static_cast<std::size_t>(d)] == Boundary::periodic) out.push_back(d);
    return out;
}

std::vector<std::vector<double>> Lattice::momenta() const {
    std::vector<int> const dirs = periodic_directions();
    if (dirs.empty()) return {std::vector<double>{}};
    // Cartesian product, last direction fastest (itertools.product).
    std::vector<std::vector<double>> rows;
    std::vector<int> j(dirs.size(), 0);
    while (true) {
        std::vector<double> row(dirs.size());
        for (std::size_t k = 0; k < dirs.size(); ++k) {
            int const L = extent_[static_cast<std::size_t>(dirs[k])];
            row[k] = 2.0 * M_PI * j[k] / L;
        }
        rows.push_back(std::move(row));
        int k = static_cast<int>(dirs.size()) - 1;
        for (; k >= 0; --k) {
            int const L = extent_[static_cast<std::size_t>(dirs[static_cast<std::size_t>(k)])];
            if (++j[static_cast<std::size_t>(k)] < L) break;
            j[static_cast<std::size_t>(k)] = 0;
        }
        if (k < 0) break;
    }
    return rows;
}

std::vector<std::vector<std::size_t>> Lattice::translation_generators() const {
    std::vector<std::vector<std::size_t>> result;
    for (int d : periodic_directions()) {
        std::size_t stride = static_cast<std::size_t>(unit_cell_.num_vertices());
        for (std::size_t axis = static_cast<std::size_t>(d) + 1; axis < extent_.size(); ++axis)
            stride *= static_cast<std::size_t>(extent_[axis]);
        auto const period = static_cast<std::size_t>(extent_[static_cast<std::size_t>(d)]);
        std::vector<std::size_t> permutation(sites_.size());
        for (std::size_t site = 0; site < sites_.size(); ++site)
            permutation[site] = (site / stride) % period + 1 == period
                ? site - (period - 1) * stride : site + stride;
        result.push_back(std::move(permutation));
    }
    return result;
}

int Lattice::num_momenta() const {
    std::vector<int> const dirs = periodic_directions();
    if (dirs.empty()) return 1;
    int total = 1;
    for (int d : dirs) total *= extent_[static_cast<std::size_t>(d)];
    return total;
}

std::vector<models::ModelEdge> Lattice::edges() const {
    std::vector<models::ModelEdge> out;
    out.reserve(bonds_.size());
    for (Bond const& b : bonds_)
        out.push_back(models::ModelEdge{b.source, b.target, b.bond_type});
    return out;
}

std::vector<int> Lattice::site_types() const {
    std::vector<int> out;
    out.reserve(sites_.size());
    for (Site const& s : sites_) out.push_back(s.site_type);
    return out;
}

// -- named factories -----------------------------------------------------------

namespace {
inline int or_default(int v, int fallback) { return v > 0 ? v : fallback; }
inline double or_default(double v, double fallback) { return v > 0.0 ? v : fallback; }
} // namespace

Lattice chain(int L, double a, Boundary bc) {
    UnitCell uc(1);
    uc.add_vertex({0.0});
    uc.add_edge(0, 0, 0, {}, {1});
    // lattices.toml ships two LATTICEGRAPHs hardcoding their BOUNDARY;
    // select the graph name to match `bc` (bc1 is not consulted).
    std::string const name = bc == Boundary::open ? "open chain lattice" : "chain lattice";
    return Lattice::build(uc, {{a}}, {L}, uniform_bc(bc, 1), name);
}

Lattice nnn_chain(int L, double a, Boundary bc) {
    UnitCell uc(1);
    uc.add_vertex({0.0});
    uc.add_edge(0, 0, 0, {}, {1});  // NN
    uc.add_edge(0, 0, 1, {}, {2});  // NNN
    std::string const name =
        bc == Boundary::periodic ? "nnn chain lattice" : "nnn open chain lattice";
    return Lattice::build(uc, {{a}}, {L}, uniform_bc(bc, 1), name);
}

Lattice square(int L, int W, double a, Boundary bc) {
    W = or_default(W, L);
    UnitCell uc(2);
    uc.add_vertex({0.0, 0.0});
    uc.add_edge(0, 0, 0, {}, {1, 0});
    uc.add_edge(0, 0, 0, {}, {0, 1});
    return Lattice::build(uc, {{a, 0}, {0, a}}, {L, W}, uniform_bc(bc, 2),
                          "square lattice");
}

Lattice frustrated_square(int L, int W, double a, Boundary bc) {
    W = or_default(W, L);
    UnitCell uc(2);
    uc.add_vertex({0.0, 0.0});
    uc.add_edge(0, 0, 0, {}, {1, 0});
    uc.add_edge(0, 0, 0, {}, {0, 1});
    uc.add_edge(0, 0, 1, {}, {1, 1});
    uc.add_edge(0, 0, 1, {}, {1, -1});
    return Lattice::build(uc, {{a, 0}, {0, a}}, {L, W}, uniform_bc(bc, 2),
                          "frustrated square lattice");
}

Lattice triangular(int L, int W, double a, Boundary bc) {
    W = or_default(W, L);
    UnitCell uc(2);
    uc.add_vertex({0.0, 0.0});
    uc.add_edge(0, 0, 0, {}, {1, 0});
    uc.add_edge(0, 0, 0, {}, {0, 1});
    uc.add_edge(0, 0, 0, {}, {-1, 1});
    return Lattice::build(uc, {{a, 0.0}, {a / 2, a * std::sqrt(3.0) / 2}}, {L, W},
                          uniform_bc(bc, 2), "triangular lattice");
}

Lattice honeycomb(int L, int W, double a, Boundary bc) {
    W = or_default(W, L);
    UnitCell uc(2);
    uc.add_vertex({0.0, 0.0}, 0);
    uc.add_vertex({-1.0 / 3.0, 2.0 / 3.0}, 1);
    uc.add_edge(0, 1, 0);                    // intra-cell
    uc.add_edge(1, 0, 1, {}, {0, 1});        // inter-cell y
    uc.add_edge(0, 1, 2, {}, {1, -1});       // inter-cell diagonal
    return Lattice::build(uc, {{a, 0.0}, {a / 2, a * std::sqrt(3.0) / 2}}, {L, W},
                          uniform_bc(bc, 2), "honeycomb lattice");
}

Lattice kagome(int L, int W, double a, Boundary bc) {
    W = or_default(W, L);
    UnitCell uc(2);
    uc.add_vertex({0.0, 0.0}, 0);
    uc.add_vertex({0.5, 0.0}, 1);
    uc.add_vertex({0.0, 0.5}, 2);
    uc.add_edge(0, 1, 0);
    uc.add_edge(0, 2, 0);
    uc.add_edge(1, 2, 0);
    uc.add_edge(0, 1, 0, {}, {-1, 0});
    uc.add_edge(1, 2, 0, {}, {1, -1});
    uc.add_edge(0, 2, 0, {}, {0, -1});
    return Lattice::build(uc, {{a, 0.0}, {a / 2, a * std::sqrt(3.0) / 2}}, {L, W},
                          uniform_bc(bc, 2), "kagome lattice");
}

Lattice simple_cubic(int L, int W, int H, double a, Boundary bc) {
    W = or_default(W, L);
    H = or_default(H, W);
    UnitCell uc(3);
    uc.add_vertex({0.0, 0.0, 0.0});
    uc.add_edge(0, 0, 0, {}, {1, 0, 0});
    uc.add_edge(0, 0, 0, {}, {0, 1, 0});
    uc.add_edge(0, 0, 0, {}, {0, 0, 1});
    return Lattice::build(uc, {{a, 0, 0}, {0, a, 0}, {0, 0, a}}, {L, W, H},
                          uniform_bc(bc, 3), "simple cubic lattice");
}

Lattice ladder(int L, double a, Boundary bc) {
    // An L x 2 square tiling so the graph matches the ALPS "ladder"
    // LATTICEGRAPH exactly (legs = x = type 0, rungs = y = type 1; the
    // rung direction is always open; W emits as 2 for the legacy bridge).
    UnitCell uc(2);
    uc.add_vertex({0.0, 0.0});
    uc.add_edge(0, 0, 0, {}, {1, 0});  // legs
    uc.add_edge(0, 0, 1, {}, {0, 1});  // rungs
    return Lattice::build(uc, {{a, 0}, {0, a}}, {L, 2}, {bc, Boundary::open},
                          "ladder");
}

Lattice rectangular(int L, int W, double a, double b, Boundary bc) {
    W = or_default(W, L);
    b = or_default(b, a);
    UnitCell uc(2);
    uc.add_vertex({0.0, 0.0});
    uc.add_edge(0, 0, 0, {}, {1, 0});
    uc.add_edge(0, 0, 0, {}, {0, 1});
    return Lattice::build(uc, {{a, 0}, {0, b}}, {L, W}, uniform_bc(bc, 2),
                          "rectangular lattice");
}

Lattice face_centered_cubic(int L, int W, int H, double a, Boundary bc) {
    W = or_default(W, L);
    H = or_default(H, W);
    UnitCell uc(3);
    uc.add_vertex({0.0, 0.0, 0.0});
    uc.add_edge(0, 0, 0, {}, {1, 0, 0});
    uc.add_edge(0, 0, 0, {}, {0, 1, 0});
    uc.add_edge(0, 0, 0, {}, {0, 0, 1});
    uc.add_edge(0, 0, 0, {}, {1, -1, 0});
    uc.add_edge(0, 0, 0, {}, {1, 0, -1});
    uc.add_edge(0, 0, 0, {}, {0, 1, -1});
    return Lattice::build(uc,
                          {{a / 2, a / 2, 0}, {a / 2, 0, a / 2}, {0, a / 2, a / 2}},
                          {L, W, H}, uniform_bc(bc, 3), "face-centered cubic lattice");
}

Lattice body_centered_cubic(int L, int W, int H, double a, Boundary bc) {
    W = or_default(W, L);
    H = or_default(H, W);
    UnitCell uc(3);
    uc.add_vertex({0.0, 0.0, 0.0});
    uc.add_edge(0, 0, 0, {}, {1, 0, 0});
    uc.add_edge(0, 0, 0, {}, {0, 1, 0});
    uc.add_edge(0, 0, 0, {}, {0, 0, 1});
    uc.add_edge(0, 0, 0, {}, {1, 1, 1});
    return Lattice::build(uc,
                          {{a / 2, a / 2, -a / 2}, {a / 2, -a / 2, a / 2},
                           {-a / 2, a / 2, a / 2}},
                          {L, W, H}, uniform_bc(bc, 3), "body-centered cubic lattice");
}

Lattice tetragonal(int L, int W, int H, double a, double c, Boundary bc) {
    W = or_default(W, L);
    H = or_default(H, W);
    UnitCell uc(3);
    uc.add_vertex({0.0, 0.0, 0.0});
    uc.add_edge(0, 0, 0, {}, {1, 0, 0});
    uc.add_edge(0, 0, 0, {}, {0, 1, 0});
    uc.add_edge(0, 0, 0, {}, {0, 0, 1});
    return Lattice::build(uc, {{a, 0, 0}, {0, a, 0}, {0, 0, c}}, {L, W, H},
                          uniform_bc(bc, 3), "tetragonal lattice");
}

Lattice orthorhombic(int L, int W, int H, double a, double b, double c, Boundary bc) {
    W = or_default(W, L);
    H = or_default(H, W);
    UnitCell uc(3);
    uc.add_vertex({0.0, 0.0, 0.0});
    uc.add_edge(0, 0, 0, {}, {1, 0, 0});
    uc.add_edge(0, 0, 0, {}, {0, 1, 0});
    uc.add_edge(0, 0, 0, {}, {0, 0, 1});
    return Lattice::build(uc, {{a, 0, 0}, {0, b, 0}, {0, 0, c}}, {L, W, H},
                          uniform_bc(bc, 3), "orthorhombic lattice");
}

Lattice hexagonal(int L, int W, int H, double a, double c, Boundary bc) {
    W = or_default(W, L);
    H = or_default(H, W);
    UnitCell uc(3);
    uc.add_vertex({0.0, 0.0, 0.0});
    uc.add_edge(0, 0, 0, {}, {1, 0, 0});
    uc.add_edge(0, 0, 0, {}, {0, 1, 0});
    uc.add_edge(0, 0, 0, {}, {0, 0, 1});
    return Lattice::build(uc,
                          {{a, 0, 0}, {a / 2, a * std::sqrt(3.0) / 2, 0}, {0, 0, c}},
                          {L, W, H}, uniform_bc(bc, 3), "hexagonal lattice");
}

} // namespace lattice
} // namespace alps
