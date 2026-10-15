// Copyright (C) 1994-2010 by Matthias Troyer, Andreas Honecker and Ryo IGARASHI
// Copyright (C) 2026 by the ALPS collaboration
// Part of the ALPS Project — see LICENSE.txt for full license text.
// SPDX-License-Identifier: MIT
//
// Ported from skilledwolf/alpsx (private repository), branch
// agent/v3-execution-and-api-audit at f8ee6e98b, file
// src/model/include/alps/model/catalog.hpp, itself extracted from this
// tree's hamiltonian_matrix.hpp (local_matrix and the fermionic bond-term
// rule, hence the 1994-2010 notice). Adaptation: legacy <alps/model.h> and
// <alps/lattice.h> headers, alps::multi_array instead of alpsx's mdarray,
// <alps/models/compiled.hpp> operators; conserved_quantum_numbers, sector
// and local_charge follow the charge derivation of alpsx's
// src/solvers/diag/src/runtime/catalog_problem.cpp.
#ifndef ALPS_MODEL_CATALOG_HPP
#define ALPS_MODEL_CATALOG_HPP

#include <alps/model.h>
#include <alps/lattice.h>
#include <alps/models/compiled.hpp>
#include <alps/multi_array.hpp>
#include <alps/numeric/is_nonzero.hpp>
#include <alps/type_traits/is_symbolic.hpp>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace alps::models {

/// Catalog expression adapter: the XML model and lattice libraries resolved
/// through the legacy graph_helper/model_helper pair and compiled to detached
/// numeric operators (OUT-row/IN-column local factors) for ALPS::ed and any
/// other consumer of CompiledOperator. It is independent of basis enumeration,
/// momentum projection and matrix storage; parity strings require
/// local_parity() when applying the returned operator. It retains the
/// catalog's cutoff, expression order and fermionic tensor conventions. This
/// is an evaluator, not a Hermiticity/conservation validator: the solver
/// applies those restrictions to the compiled operator. Lazy local-basis
/// initialization and parameter updates require external synchronization if
/// one instance is shared between threads.
template <class Scalar = std::complex<double>, class G = typename graph_helper<>::graph_type>
class CatalogModel {
public:
  using value_type = Scalar;
  using site_descriptor = typename graph_helper<G>::site_descriptor;
  using bond_descriptor = typename graph_helper<G>::bond_descriptor;

  explicit CatalogModel(Parameters const& parameters)
    : parms(parameters), graph_(parameters),
      model_(graph_, parameters, is_symbolic<Scalar>::type::value) {}
  // graph_helper retains internal property-map references. Transfer the
  // compiled values, rather than copying or relocating the evaluator.
  CatalogModel(CatalogModel const&) = delete;
  CatalogModel& operator=(CatalogModel const&) = delete;
  CatalogModel(CatalogModel&&) = delete;
  CatalogModel& operator=(CatalogModel&&) = delete;

  /// Merge evaluation/sector parameters, invalidating the local basis. Graph
  /// topology and resolved catalog definitions remain those of construction;
  /// construct another model to change the lattice or model definition.
  void set_parameters(Parameters const& parameters) { parms << parameters; built_basis_ = false; }
  Parameters const& parameters() const { return parms; }
  graph_helper<G> const& graph() const { return graph_; }
  model_helper<> const& model() const { return model_; }

  /// Local states only: no product states, Bloch orbits or matrix are allocated.
  basis_states_descriptor<short> const& basis() const {
    if (!built_basis_) {
      typename model_helper<>::basis_descriptor_type descriptor = model_.basis();
      descriptor.set_parameters(parms);
      basis_ = basis_states_descriptor<short>(descriptor, graph_.graph());
      built_basis_ = true;
    }
    return basis_;
  }

  /// Does not reseed the expression RNG; callers control their Disorder context.
  /// Site, bond, global, named and explicit two-site tensor operators share this
  /// compiler. Optional locations select a site, bond, or ordered pair of sites.
  template <class Op, class... Location>
  CompiledOperator compile_operator(Op const& op, Location... location) const {
    static_assert(std::is_same_v<Scalar, double> || std::is_same_v<Scalar, std::complex<double>>,
                  "catalog numeric compilation requires double or complex<double>");
    CompiledOperator result;
    for (auto const& local : basis()) result.site_dims.push_back(static_cast<int>(local.size()));
    append_operator(result, op, location...);
    return result;
  }
  CompiledOperator compile_hamiltonian() const { return compile_operator(model_.model()); }
  /// One +1/-1 diagonal per site: -1 on local states with odd fermion number.
  std::vector<std::vector<double>> local_parity() const;

  /// Conserved charges, sorted: CONSERVED_QUANTUMNUMBERS plus every constraint
  /// the parameters evaluate (Sz_total, N_total, ...). Each must be a
  /// quantum number constrained by the model's basis.
  std::vector<std::string> conserved_quantum_numbers() const;
  /// The evaluated constraints as totals keyed by quantum number.
  std::map<std::string, double> sector() const;
  /// Per site, per local state value of one quantum number, in basis() order.
  std::vector<std::vector<double>> local_charge(std::string const& name) const;

  multi_array<value_type,2> local_matrix(SiteOperator const&, site_descriptor) const;
  multi_array<std::pair<value_type,bool>,4> local_matrix(BondOperator const&, bond_descriptor const&) const;
  multi_array<std::pair<value_type,bool>,4> local_matrix(BondOperator const&, site_descriptor const&, site_descriptor const&) const;

private:
  void append_operator(CompiledOperator&, SiteOperator const&, site_descriptor) const;
  void append_operator(CompiledOperator&, BondOperator const&, bond_descriptor) const;
  void append_operator(CompiledOperator&, BondOperator const&, site_descriptor, site_descriptor) const;
  void append_operator(CompiledOperator&, multi_array<std::pair<value_type,bool>,4> const&, site_descriptor, site_descriptor) const;
  void append_operator(CompiledOperator&, SiteOperator const&) const;
  void append_operator(CompiledOperator&, BondOperator const&) const;
  void append_operator(CompiledOperator&, GlobalOperator const&) const;
  void append_operator(CompiledOperator&, std::string const&) const;
  void append_operator(CompiledOperator&, std::string const&, site_descriptor) const;
  void append_operator(CompiledOperator&, std::string const&, bond_descriptor) const;

  Parameters parms;
  graph_helper<G> graph_;
  model_helper<> model_;
  mutable bool built_basis_ = false;
  mutable basis_states_descriptor<short> basis_;
};

template <class Scalar, class G>
multi_array<typename CatalogModel<Scalar,G>::value_type,2> CatalogModel<Scalar,G>::local_matrix(const SiteOperator& op, site_descriptor s) const
{
  Parameters p(parms);
  if (graph_.inhomogeneous_sites()) {
    throw_if_xyz_defined(parms,s); // check whether x, y, or z is set
    p << graph_.coordinate_as_parameter(s); // set x, y and z
  }
  return get_matrix(typename CatalogModel<Scalar,G>::value_type(),op,model_.model().basis().site_basis(graph_.site_type(s)),p);
}

template <class Scalar, class G>
multi_array<std::pair<typename CatalogModel<Scalar,G>::value_type,bool>,4> CatalogModel<Scalar,G>::local_matrix(const BondOperator& op, const bond_descriptor& b) const
{
  unsigned int stype1 = graph_.site_type(graph_.source(b));
  unsigned int stype2 = graph_.site_type(graph_.target(b));
  Parameters p(parms);
  if (graph_.inhomogeneous_bonds()) {
    throw_if_xyz_defined(parms,b); // check whether x, y, or z is set
    p << graph_.coordinate_as_parameter(b); // set x, y and z
  }
  return get_fermionic_matrix(typename CatalogModel<Scalar,G>::value_type(),op,model_.model().basis().site_basis(stype1),
                                      model_.model().basis().site_basis(stype2),p);
}

template <class Scalar, class G>
multi_array<std::pair<typename CatalogModel<Scalar,G>::value_type,bool>,4> CatalogModel<Scalar,G>::local_matrix(const BondOperator& op, const site_descriptor& s1, const site_descriptor& s2) const
{
  unsigned int stype1 = graph_.site_type(s1);
  unsigned int stype2 = graph_.site_type(s2);
  return get_fermionic_matrix(typename CatalogModel<Scalar,G>::value_type(),op,model_.model().basis().site_basis(stype1),
                                      model_.model().basis().site_basis(stype2),parms);
}

template <class Scalar, class G>
void CatalogModel<Scalar,G>::append_operator(CompiledOperator& result,
                                         const SiteOperator& op, site_descriptor site) const
{
  auto const mat = local_matrix(op, site);
  // Local tensors are input-first. Only the catalog's established cutoff is
  // applied; generic compiled observables retain arbitrary complex values.
  std::size_t const dim = mat.shape()[0];
  LocalFactor factor{static_cast<int>(site), std::vector<std::complex<double>>(dim * dim)};
  bool nonzero = false;
  for (std::size_t input = 0; input < dim; ++input)
    for (std::size_t output = 0; output < dim; ++output)
      if (numeric::is_nonzero(mat[input][output])) {
        factor.matrix[output * dim + input] = mat[input][output];
        nonzero = true;
      }
  if (nonzero) result.terms.push_back({1., {std::move(factor)}, {}});
}

template <class Scalar, class G>
void CatalogModel<Scalar,G>::append_operator(CompiledOperator& result,
                                         const BondOperator& op, bond_descriptor bond) const
{
  // The bond descriptor injects midpoint coordinates through local_matrix.
  append_operator(result, local_matrix(op, bond), graph_.source(bond), graph_.target(bond));
}

template <class Scalar, class G>
void CatalogModel<Scalar,G>::append_operator(CompiledOperator& result,
                                         const BondOperator& op, site_descriptor s1, site_descriptor s2) const
{
  // Arbitrary endpoint pairs retain the original parameter-only evaluation.
  append_operator(result, local_matrix(op, s1, s2), s1, s2);
}

template <class Scalar, class G>
void CatalogModel<Scalar,G>::append_operator(CompiledOperator& result,
                                         const multi_array<std::pair<value_type,bool>,4>& mat,
                                         site_descriptor s1, site_descriptor s2) const
{
  std::size_t const d1 = mat.shape()[0], d2 = mat.shape()[1];
  if (s1 == s2) {
    // The existing catalog rule writes endpoint two last. Keep a separate
    // term for every first-endpoint output to preserve accumulation order.
    for (std::size_t j1 = 0; j1 < d1; ++j1) {
      LocalFactor factor{static_cast<int>(s1), std::vector<std::complex<double>>(d1 * d1)};
      bool nonzero = false;
      for (std::size_t input = 0; input < d1; ++input)
        for (std::size_t output = 0; output < d2; ++output) {
          auto const& entry = mat[input][input][j1][output];
          if (numeric::is_nonzero(entry.first)) {
            factor.matrix[output * d1 + input] = entry.second ? -entry.first : entry.first;
            nonzero = true;
          }
        }
      if (nonzero) result.terms.push_back({1., {std::move(factor)}, {}});
    }
    return;
  }
  // The legacy sign rule: a fermionic tensor picks up (-1) when the second
  // endpoint follows the first, times the parity of the source state on every
  // site from min(s1,s2) up to but excluding max(s1,s2). The same rule is
  // expressed as a coefficient sign plus a Jordan-Wigner string that includes
  // the lower endpoint, whose parity ALPS::ed reads from the input state.
  for (std::size_t i1 = 0; i1 < d1; ++i1)
    for (std::size_t j1 = 0; j1 < d1; ++j1)
      for (bool fermionic : {false, true}) {
        LocalFactor second{static_cast<int>(s2), std::vector<std::complex<double>>(d2 * d2)};
        bool nonzero = false;
        for (std::size_t i2 = 0; i2 < d2; ++i2)
          for (std::size_t j2 = 0; j2 < d2; ++j2) {
            auto const& entry = mat[i1][i2][j1][j2];
            if (entry.second == fermionic && numeric::is_nonzero(entry.first)) {
              second.matrix[j2 * d2 + i2] = entry.first;
              nonzero = true;
            }
          }
        if (!nonzero) continue;
        LocalFactor first{static_cast<int>(s1), std::vector<std::complex<double>>(d1 * d1)};
        first.matrix[j1 * d1 + i1] = 1.;
        LocalTerm term{fermionic && s2 >= s1 ? -1. : 1., {std::move(first), std::move(second)}, {}};
        if (fermionic)
          for (auto site = std::min(s1, s2); site < std::max(s1, s2); ++site)
            term.jordan_wigner.push_back(static_cast<int>(site));
        result.terms.push_back(std::move(term));
      }
}

template <class Scalar, class G>
void CatalogModel<Scalar,G>::append_operator(CompiledOperator& result, const SiteOperator& op) const
{
  for (auto it = graph_.sites().first; it != graph_.sites().second; ++it)
    append_operator(result, op, *it);
}

template <class Scalar, class G>
void CatalogModel<Scalar,G>::append_operator(CompiledOperator& result, const BondOperator& op) const
{
  for (auto it = graph_.bonds().first; it != graph_.bonds().second; ++it)
    append_operator(result, op, *it);
}

template <class Scalar, class G>
void CatalogModel<Scalar,G>::append_operator(CompiledOperator& result, const GlobalOperator& op) const
{
  // Preserve expression evaluation order, including random calls.
  for (auto it = graph_.sites().first; it != graph_.sites().second; ++it)
    append_operator(result, op.site_term(graph_.site_type(*it)), *it);
  for (auto it = graph_.bonds().first; it != graph_.bonds().second; ++it)
    append_operator(result, op.bond_term(graph_.bond_type(*it)), *it);
}

template <class Scalar, class G>
void CatalogModel<Scalar,G>::append_operator(CompiledOperator& result,
                                         const std::string& name, site_descriptor site) const
{
  if (model_.has_site_operator(name)) append_operator(result, model_.get_site_operator(name), site);
  else append_operator(result, SiteOperator(name+"(i)","i"), site);
}

template <class Scalar, class G>
void CatalogModel<Scalar,G>::append_operator(CompiledOperator& result,
                                         const std::string& name, bond_descriptor bond) const
{
  append_operator(result, model_.get_bond_operator(name), bond);
}

template <class Scalar, class G>
void CatalogModel<Scalar,G>::append_operator(CompiledOperator& result, const std::string& name) const
{
  if (model_.has_site_operator(name)) append_operator(result, model_.get_site_operator(name));
  else if (model_.has_bond_operator(name)) append_operator(result, model_.get_bond_operator(name));
  else if (model_.has_global_operator(name)) append_operator(result, model_.get_global_operator(name));
  else append_operator(result, SiteOperator(name+"(i)","i"));
}

template <class Scalar, class G>
std::vector<std::vector<double>> CatalogModel<Scalar,G>::local_parity() const
{
  auto const& local_basis = basis();
  std::vector<std::vector<double>> parity(local_basis.size());
  for (std::size_t site = 0; site < local_basis.size(); ++site)
    for (std::size_t state = 0; state < local_basis[site].size(); ++state)
      parity[site].push_back(is_fermionic(model_.model().basis().site_basis(graph_.site_type(site)),
                                         local_basis[site][state]) ? -1. : 1.);
  return parity;
}

template <class Scalar, class G>
std::vector<std::string> CatalogModel<Scalar,G>::conserved_quantum_numbers() const
{
  auto const& descriptor = basis().get_basis();
  std::set<std::string> available, declared;
  for (auto const& constraint : descriptor.all_constraints()) available.insert(constraint.first);
  for (auto const& constraint : descriptor.constraints()) declared.insert(constraint.first);
  std::string const list = parms.value_or_default("CONSERVED_QUANTUMNUMBERS", "");
  std::string const separators = " ,\t\r\n";
  for (std::size_t begin = list.find_first_not_of(separators); begin != std::string::npos;) {
    auto const end = list.find_first_of(separators, begin);
    declared.insert(list.substr(begin, end == std::string::npos ? end : end - begin));
    begin = list.find_first_not_of(separators, end);
  }
  for (auto const& name : declared)
    if (!available.count(name))
      throw std::invalid_argument("catalog model: unknown conserved quantum number " + name);
  return std::vector<std::string>(declared.begin(), declared.end());
}

template <class Scalar, class G>
std::map<std::string, double> CatalogModel<Scalar,G>::sector() const
{
  std::map<std::string, double> result;
  for (auto const& constraint : basis().get_basis().constraints())
    result.emplace(constraint.first, to_double(constraint.second));
  return result;
}

template <class Scalar, class G>
std::vector<std::vector<double>> CatalogModel<Scalar,G>::local_charge(std::string const& name) const
{
  auto const& local_basis = basis();
  std::vector<std::vector<double>> rows(local_basis.size());
  for (std::size_t site = 0; site < local_basis.size(); ++site) {
    auto const& descriptor = local_basis.get_site_basis(static_cast<int>(site));
    auto const index = get_quantumnumber_index(name, descriptor);
    if (index >= descriptor.size())
      throw std::invalid_argument("catalog model: site " + std::to_string(site) + " has no quantum number " + name);
    for (auto const& state : local_basis[site]) rows[site].push_back(to_double(get_quantumnumber(state, index)));
  }
  return rows;
}

} // namespace alps::models

#endif // ALPS_MODEL_CATALOG_HPP
