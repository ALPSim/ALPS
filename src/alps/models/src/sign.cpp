// Copyright (C) 2003-2009 by Matthias Troyer <troyer@itp.phys.ethz.ch>,
//                            Synge Todo <wistaria@comp-phys.org>
// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Adapted from this tree's src/alps/model/include/alps/model/sign.h
// (has_sign_problem, detail::BondSiteSignVisitor, is_frustrated). Changes:
// the matrices are the summed one- and two-site terms of a CompiledOperator,
// per site and per site pair instead of per site type and bond type of a
// catalog model, so the graph has one edge per interacting pair and no
// inhomogeneous-bond fallback is needed; zero-sign pairs are left out of the
// graph instead of filtered; inputs outside the rule throw; C++17.
#include <alps/models/sign.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/undirected_dfs.hpp>
#include <boost/property_map/property_map.hpp>

#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace alps {
namespace models {

namespace {

// +1 for a negative element (a positive weight of -H), -1 for a positive one,
// 0 within the legacy cutoff; false once two nonzero elements disagree.
class MatrixSign {
public:
    bool add(double value) {
        int const next = value < -1.0e-10 ? 1 : (value > 1.0e-10 ? -1 : 0);
        if (!sign_) sign_ = next;
        else if (next && sign_ != next) return false;
        return true;
    }
    int value() const { return sign_; }

private:
    int sign_ = 0;
};

// The legacy visitor, unchanged: the DFS assigns relative sublattice signs
// and flags a tree edge whose bond and site signs disagree, or a back edge
// that closes an inconsistent cycle.
template <class Graph, class PropertyMap, class BondPropertyMap, class SitePropertyMap>
class BondSiteSignVisitor : public boost::dfs_visitor<> {
public:
    using vertex_descriptor = typename boost::graph_traits<Graph>::vertex_descriptor;
    using edge_descriptor = typename boost::graph_traits<Graph>::edge_descriptor;

    BondSiteSignVisitor(PropertyMap map, bool* check, BondPropertyMap bond_sign, SitePropertyMap site_sign)
        : map_(map), check_(check), bond_sign_(bond_sign), site_sign_(site_sign) { *check_ = false; }

    void initialize_vertex(vertex_descriptor s, Graph const&) { map_[s] = 0; }
    void start_vertex(vertex_descriptor s, Graph const&) { map_[s] = 1; }
    void tree_edge(edge_descriptor e, Graph const& g) {
        if (map_[boost::target(e, g)] == 0) {
            if (bond_sign_[e] > 0) map_[boost::target(e, g)] = map_[boost::source(e, g)];
            else map_[boost::target(e, g)] = -map_[boost::source(e, g)];
        } else {
            if (bond_sign_[e] > 0) map_[boost::source(e, g)] = map_[boost::target(e, g)];
            else map_[boost::source(e, g)] = -map_[boost::target(e, g)];
        }
        if (bond_sign_[e] * site_sign_[boost::source(e, g)] * site_sign_[boost::target(e, g)] < 0) *check_ = true;
    }
    void back_edge(edge_descriptor e, Graph const& g) {
        if (bond_sign_[e] * map_[boost::source(e, g)] * map_[boost::target(e, g)] < 0
            || site_sign_[boost::source(e, g)] * site_sign_[boost::target(e, g)] * map_[boost::source(e, g)]
                       * map_[boost::target(e, g)] < 0)
            *check_ = true;
    }

private:
    PropertyMap map_;
    bool* check_;
    BondPropertyMap bond_sign_;
    SitePropertyMap site_sign_;
};

double real_part(Complex value) {
    if (!(std::abs(value.imag()) <= 1.0e-10))
        throw std::invalid_argument("has_sign_problem needs real matrix elements");
    return value.real();
}

} // namespace

bool has_sign_problem(CompiledOperator const& op) {
    auto const num_sites = op.site_dims.size();
    auto dim = [&](int site) { return static_cast<std::size_t>(op.site_dims[static_cast<std::size_t>(site)]); };
    // Summed before any sign is read: XX + YY cancels the double flips only in the sum.
    std::vector<std::vector<Complex>> site_matrix(num_sites);
    std::map<std::pair<int, int>, std::vector<Complex>> pair_matrix;
    for (auto const& term : op.terms) {
        if (!term.jordan_wigner.empty())
            throw std::invalid_argument("has_sign_problem does not cover Jordan-Wigner strings (fermion hopping)");
        if (term.factors.size() > 2)
            throw std::invalid_argument("has_sign_problem covers one- and two-site terms only");
        for (auto const& factor : term.factors) {
            if (factor.site < 0 || std::size_t(factor.site) >= num_sites)
                throw std::invalid_argument("has_sign_problem: factor site outside the operator");
            if (factor.matrix.size() != dim(factor.site) * dim(factor.site))
                throw std::invalid_argument("has_sign_problem: factor matrix does not match its site");
        }
        if (term.factors.size() == 1) {
            auto const& f = term.factors[0];
            auto& m = site_matrix[std::size_t(f.site)];
            if (m.empty()) m.assign(f.matrix.size(), Complex(0.));
            for (std::size_t k = 0; k < m.size(); ++k) m[k] += term.coefficient * f.matrix[k];
        } else if (term.factors.size() == 2) {
            auto const* a = &term.factors[0];
            auto const* b = &term.factors[1];
            if (a->site == b->site) throw std::invalid_argument("has_sign_problem: two factors on one site");
            if (a->site > b->site) std::swap(a, b);
            std::size_t const da = dim(a->site), db = dim(b->site), d = da * db;
            auto& m = pair_matrix[{a->site, b->site}];
            if (m.empty()) m.assign(d * d, Complex(0.));
            // Row (i, j) is the output state of (a, b), column (k, l) the input state.
            for (std::size_t i = 0; i < da; ++i)
                for (std::size_t k = 0; k < da; ++k) {
                    Complex const ak = term.coefficient * a->matrix[i * da + k];
                    if (ak == Complex(0.)) continue;
                    for (std::size_t j = 0; j < db; ++j)
                        for (std::size_t l = 0; l < db; ++l)
                            m[(i * db + j) * d + k * db + l] += ak * b->matrix[j * db + l];
                }
        }
    }
    // Every element must be real before any sign decides the answer.
    for (auto const& m : site_matrix)
        for (auto const& value : m) real_part(value);
    for (auto const& [pair, m] : pair_matrix)
        for (auto const& value : m) real_part(value);

    std::vector<int> site_sign(num_sites, 0);
    for (std::size_t s = 0; s < num_sites; ++s) {
        auto const& m = site_matrix[s];
        if (m.empty()) continue;
        std::size_t const d = dim(int(s));
        MatrixSign sign;
        for (std::size_t i = 0; i < d; ++i)
            for (std::size_t j = 0; j < d; ++j)
                if (i != j && !sign.add(m[i * d + j].real())) return true;  // indefinite site matrix
        site_sign[s] = sign.value();
    }

    using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, boost::no_property,
                                        boost::property<boost::edge_index_t, std::size_t>>;
    Graph graph(num_sites);
    std::vector<int> bond_sign;
    for (auto const& [pair, m] : pair_matrix) {
        std::size_t const da = dim(pair.first), db = dim(pair.second), d = da * db;
        MatrixSign sign;
        for (std::size_t i = 0; i < da; ++i)
            for (std::size_t j = 0; j < db; ++j)
                for (std::size_t k = 0; k < da; ++k)
                    for (std::size_t l = 0; l < db; ++l)
                        if (i != k && j != l && !sign.add(m[(i * db + j) * d + k * db + l].real()))
                            return true;  // indefinite bond matrix
        if (!sign.value()) continue;  // no element changes both sites: no edge, as the legacy filter
        auto const edge = boost::add_edge(std::size_t(pair.first), std::size_t(pair.second), graph).first;
        boost::put(boost::edge_index, graph, edge, bond_sign.size());
        bond_sign.push_back(sign.value());
    }

    std::vector<int> sublattice(num_sites, 0);
    auto const vertex_index = boost::get(boost::vertex_index, graph);
    auto const edge_index = boost::get(boost::edge_index, graph);
    auto const sublattice_map = boost::make_iterator_property_map(sublattice.begin(), vertex_index);
    auto const bond_map = boost::make_iterator_property_map(bond_sign.begin(), edge_index);
    auto const site_map = boost::make_iterator_property_map(site_sign.begin(), vertex_index);
    bool check = false;
    BondSiteSignVisitor<Graph, decltype(sublattice_map), decltype(bond_map), decltype(site_map)> visitor(
        sublattice_map, &check, bond_map, site_map);
    std::vector<boost::default_color_type> vertex_color(num_sites);
    std::vector<boost::default_color_type> edge_color(bond_sign.size());
    boost::undirected_dfs(graph, visitor, boost::make_iterator_property_map(vertex_color.begin(), vertex_index),
                          boost::make_iterator_property_map(edge_color.begin(), edge_index));
    return check;
}

bool has_sign_problem(Hamiltonian const& model) { return has_sign_problem(compile_hamiltonian(model).op); }

} // namespace models
} // namespace alps
