/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@comp-phys.org>,
*                            Andreas Honecker <ahoneck@uni-goettingen.de>,
*                            Ryo IGARASHI <rigarash@hosi.phys.s.u-tokyo.ac.jp>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_MODEL_HAMILTONIAN_MATRIX_HPP
#define ALPS_MODEL_HAMILTONIAN_MATRIX_HPP

#include <alps/model.h>
#include <alps/lattice.h>

#include <alps/numeric/numeric_cast.hpp>
#include <alps/numeric/matrix_as_vector.hpp>
#include <alps/numeric/conj.hpp>
#include <alps/numeric/is_nonzero.hpp>
#include <alps/type_traits/is_symbolic.hpp>
#include <alps/multi_array.hpp>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>

namespace alps {

namespace numeric {
    template <typename T, typename MemoryBlock>
    class matrix;
} // end namespace numeric

namespace detail {

    // The alps::numeric::matrix and ublas matrices behave differently on construction.
    // - ublas matrices are uninitalized -> call clear() afterwards
    // - alps::numeric::matrix is initialized after construction
    template <typename Matrix>
    Matrix initialized_matrix(Matrix const&, std::size_t const rows, std::size_t const cols)
    {
        Matrix m(rows,cols);
        m.clear();

        // Boost bug workaround
        // force at least one entry per row for sparse ublas matrices
        assert(rows == cols);
        for (std::size_t i=0;i<rows;++i) {
          m(i,i) +=1.;
          m(i,i) -=1.;
        }
        return m;
    }

    template <typename T, typename MemoryBlock>
    alps::numeric::matrix<T,MemoryBlock> initialized_matrix(alps::numeric::matrix<T,MemoryBlock> const&, std::size_t const rows, std::size_t const cols)
    {
        // the dense storage has rows*cols elements; a wrapped product would
        // allocate too little
        if (cols && rows > std::numeric_limits<std::size_t>::max() / cols)
          boost::throw_exception(std::runtime_error("Hilbert space dimension "
            + boost::lexical_cast<std::string>(rows) + ": a dense " + boost::lexical_cast<std::string>(rows)
            + " x " + boost::lexical_cast<std::string>(cols) + " matrix has more elements than a std::size_t can count."));
        return alps::numeric::matrix<T,MemoryBlock>(rows,cols);
    }
} // end namespace detail

template <class M, class G = typename graph_helper<>::graph_type>
class hamiltonian_matrix
{
public:
  typedef M matrix_type;
  typedef typename matrix_type::value_type value_type;
  typedef typename graph_helper<>::graph_type graph_type;
  typedef typename graph_helper<>::site_descriptor site_descriptor;
  typedef typename graph_helper<>::bond_descriptor bond_descriptor;
  typedef typename graph_helper<>::site_iterator site_iterator;
  typedef typename graph_helper<>::bond_iterator bond_iterator;
  typedef typename graph_helper<>::vector_type vector_type;
  typedef typename model_helper<>::basis_descriptor_type basis_descriptor_type;
  typedef basis_states<short> basis_states_type;
  typedef bloch_basis_states<short> bloch_basis_states_type;
  typedef basis_states_type::value_type state_type;
  
  hamiltonian_matrix (Parameters const& parms);
  void set_parameters(Parameters const& p) { parms << p ; built_basis_=false; built_matrix_=false;}
  basis_states_type& states_vector() { if (!built_basis_) build_basis(); return states; }
  const basis_states_type& states_vector() const {if (!built_basis_) build_basis(); return states; }
  bloch_basis_states_type& bloch_states_vector() { if (!built_basis_) build_basis(); return bloch_states; }
  const bloch_basis_states_type& bloch_states_vector() const { if (!built_basis_) build_basis(); return bloch_states; }
  matrix_type& matrix() {if (!built_matrix_) build(); return matrix_;}
  const matrix_type& matrix() const {if (!built_matrix_) build(); return matrix_;}
  std::size_t dimension() const { if (!built_basis_) build_basis(); return uses_translation_invariance() ? bloch_states.size() : states.size();}
  void dostep();
  void print_basis(std::ostream& os) const
  {
    if (uses_translation_invariance())
      os << bloch_states_vector();
    else
      os << states_vector();
  }


  template <class STATES, class V, class W>
  void apply_operator(const STATES&, const SiteOperator& op, site_descriptor s, const V&, W&) const;

  template <class V, class W> 
  void apply_operator(const SiteOperator& op, site_descriptor s, const V& x, W& y) const
  {
    if (uses_translation_invariance())
      apply_operator(bloch_states,op,s,x,y);
    else
      apply_operator(states,op,s,x,y);
  }
  
  template <class V, class W> 
  void apply_operator(const BondOperator& op, bond_descriptor b, const V&, W&) const;

  template <class V, class W> 
  void apply_operator(const BondOperator& op, site_descriptor s1, site_descriptor s2, const V&, W&) const;

  template <class V, class W> 
  void apply_operator(const boost::multi_array<std::pair<value_type,bool>,4>& mat, site_descriptor s1, site_descriptor s2, const V& x, W& y) const
  {
    if (uses_translation_invariance())
      apply_operator(bloch_states,mat,s1,s2,x,y);
    else
      apply_operator(states,mat,s1,s2,x,y);
  }

  template <class STATES, class V, class W>
  void apply_operator(const STATES&, const boost::multi_array<std::pair<value_type,bool>,4>& mat, site_descriptor s1, site_descriptor s2, const V& x, W& y) const;

  // n-site term on one n-site group (an ordered tuple of sites)
  template <class V, class W>
  void apply_operator(const NSiteOperator& op, const nsite_group& p, const V& x, W& y) const
  { apply_operator(local_matrix(op,p),p,x,y); }

  template <class V, class W>
  void apply_operator(const nsite_matrix<value_type>& mat, const nsite_group& p, const V& x, W& y) const
  {
    if (uses_translation_invariance())
      apply_operator(bloch_states,mat,p,x,y);
    else
      apply_operator(states,mat,p,x,y);
  }

  template <class STATES, class V, class W>
  void apply_operator(const STATES&, const nsite_matrix<value_type>& mat, const nsite_group& p, const V& x, W& y) const;

  template <class V, class W> 
  void apply_operator(const SiteOperator& op, const V&, W&) const;
  
  template <class V, class W> 
  void apply_operator(const BondOperator& op, const V&, W&) const;

  template <class V, class W> 
  void apply_operator(const GlobalOperator& op, const V&, W&) const;

  template <class MM, class OP> 
  MM operator_matrix(const OP& op) const 
  {
    MM m(detail::initialized_matrix(MM(),dimension(),dimension()));
    add_operator_matrix(m,op);
    return m;
  }

 
  template <class MM, class OP> 
  void add_operator_matrix(MM& m, const OP& op) const 
  {
    numeric::matrix_as_vector<MM> v(m);
    apply_operator(op,v,v);
  }

 
  template <class MM, class OP, class D> 
  MM operator_matrix(const OP& op, D d) const 
  {
    MM m(detail::initialized_matrix(MM(),dimension(),dimension()));
    add_operator_matrix(m,op,d);
    return m;
  }

  template <class MM, class OP, class D> 
  void add_operator_matrix(MM& m, const OP& op, D d) const 
  {
    numeric::matrix_as_vector<MM> v(m);
    apply_operator(op,d,v,v);
  }

  template <class MM, class OP> 
  MM operator_matrix(const OP& op, site_descriptor s1, site_descriptor s2) const 
  {
    MM m(detail::initialized_matrix(MM(),dimension(),dimension()));
    add_operator_matrix(m,op,s1,s2);
    return m;
  }

  template <class MM, class OP> 
  void add_operator_matrix(MM& m, const OP& op, site_descriptor s1, site_descriptor s2) const 
  {
    numeric::matrix_as_vector<MM> v(m);
    apply_operator(op,s1,s2,v,v);
  }

  multi_array<value_type,2> local_matrix(const SiteOperator& op, site_descriptor s) const;
  multi_array<std::pair<value_type,bool>,4> local_matrix(const BondOperator& op, const bond_descriptor& b) const;
  multi_array<std::pair<value_type,bool>,4> local_matrix(const BondOperator& op, const site_descriptor& s1, const site_descriptor& s2) const;
  nsite_matrix<value_type> local_matrix(const NSiteOperator& op, const nsite_group& p) const;

  template <class V, class W> 
  void apply_operator(const std::string& name, bond_descriptor b, const V& x, W& y) const
  { apply_operator(model_.get_bond_operator(name),b,x,y); }

  template <class V, class W>
  void apply_operator(const std::string& op, site_descriptor s, const V& x, W& y) const;

  template <class V, class W> 
  void apply_operator(const std::string& name, const V& x, W& y) const;

  bool uses_translation_invariance() const { return parms.defined("TOTAL_MOMENTUM") && !static_cast<std::string>(parms["TOTAL_MOMENTUM"]).empty();}

protected:
  void build() const;
  void build_basis() const;
  void warn_if_cutoff(const nsite_matrix<value_type>& mat, const nsite_group& p) const;
  void check_sector(const state_type& state, const state_type& newstate,
                    const std::string& term, const std::vector<std::size_t>& sites) const;

  mutable basis_states_type states;
  mutable bloch_basis_states_type bloch_states;

private:
  Parameters parms;
  mutable bool built_matrix_;
  mutable bool built_basis_;
  mutable bool warned_cutoff_;
  // true while build() fills the Hamiltonian; operator_matrix() for
  // measurements leaves it false
  mutable bool building_hamiltonian_;
  mutable matrix_type matrix_;
  mutable basis_states_descriptor<short> basis_;
  graph_helper<G> graph_;
  model_helper<> model_;
  
};


template <class M, class G>
hamiltonian_matrix<M,G>::hamiltonian_matrix(Parameters const& p)
  : parms(p)
  ,  built_matrix_(false)
  ,  built_basis_(false)
  ,  warned_cutoff_(false)
  ,  building_hamiltonian_(false)
  ,  graph_(p)
  ,  model_(graph_,p,is_symbolic<value_type>::type::value,true)
{}    


template <class M, class G>
multi_array<typename hamiltonian_matrix<M,G>::value_type,2> hamiltonian_matrix<M,G>::local_matrix(const SiteOperator& op, site_descriptor s) const
{
  Parameters p(parms);
  if (graph_.inhomogeneous_sites()) {
    throw_if_xyz_defined(parms,s); // check whether x, y, or z is set
    p << graph_.coordinate_as_parameter(s); // set x, y and z
  }
  return get_matrix(typename hamiltonian_matrix<M,G>::value_type(),op,model_.model().basis().site_basis(graph_.site_type(s)),p);
}

template <class M, class G>
multi_array<std::pair<typename hamiltonian_matrix<M,G>::value_type,bool>,4> hamiltonian_matrix<M,G>::local_matrix(const BondOperator& op, const bond_descriptor& b) const
{
  unsigned int stype1 = graph_.site_type(graph_.source(b));
  unsigned int stype2 = graph_.site_type(graph_.target(b));
  Parameters p(parms);
  if (graph_.inhomogeneous_bonds()) {
    throw_if_xyz_defined(parms,b); // check whether x, y, or z is set
    p << graph_.coordinate_as_parameter(b); // set x, y and z
  }
  return get_fermionic_matrix(typename hamiltonian_matrix<M,G>::value_type(),op,model_.model().basis().site_basis(stype1),
                                      model_.model().basis().site_basis(stype2),p);  
}

template <class M, class G>
multi_array<std::pair<typename hamiltonian_matrix<M,G>::value_type,bool>,4> hamiltonian_matrix<M,G>::local_matrix(const BondOperator& op, const site_descriptor& s1, const site_descriptor& s2) const
{
  unsigned int stype1 = graph_.site_type(s1);
  unsigned int stype2 = graph_.site_type(s2);
  return get_fermionic_matrix(typename hamiltonian_matrix<M,G>::value_type(),op,model_.model().basis().site_basis(stype1),
                                      model_.model().basis().site_basis(stype2),parms);  
}

template <class M, class G>
nsite_matrix<typename hamiltonian_matrix<M,G>::value_type> hamiltonian_matrix<M,G>::local_matrix(const NSiteOperator& op, const nsite_group& p) const
{
  std::vector<SiteBasisDescriptor<short> > bases;
  for (std::size_t r=0;r<p.sites.size();++r)
    bases.push_back(model_.model().basis().site_basis(graph_.site_type(p.sites[r])));
  return op.template matrix<value_type>(bases,parms);
}

// A term that changes a quantum number whose range is capped by a parameter
// (e.g. a boson N with max="Nmax") gives results that depend on that cutoff:
// entries that would leave the range are dropped. Say so once. Fixed ranges
// (fermion 0..1) and ranges set by another quantum number (Sz in -S..S) are
// physical, not a cutoff.
template <class M, class G>
void hamiltonian_matrix<M,G>::warn_if_cutoff(const nsite_matrix<value_type>& mat, const nsite_group& p) const
{
  if (warned_cutoff_ || mat.nonconserved.empty())
    return;
  for (std::size_t r=0;r<p.sites.size();++r) {
    SiteBasisDescriptor<short> b = model_.model().basis().site_basis(graph_.site_type(p.sites[r]));
    b.set_parameters(parms);
    for (std::size_t q=0;q<b.size();++q) {
      if (!mat.nonconserved.count(b[q].name()) || can_evaluate(b[q].max_expression()))
        continue;
      bool set_by_quantumnumber=false;
      for (std::size_t o=0;o<b.size();++o)
        if (o!=q && b[q].depends_on(b[o].name()))
          set_by_quantumnumber=true;
      if (set_by_quantumnumber)
        continue;
      std::cerr << "Warning: an <NSITETERM> changes the total " << b[q].name()
                << ", whose local range is cut off at " << b[q].name() << " <= "
                << b[q].max_expression() << " = " << b[q].max() << ". Matrix elements "
                << "that would exceed the cutoff are dropped, so results depend on it; "
                << "check convergence in " << b[q].max_expression() << ".\n";
      warned_cutoff_=true;
      return;
    }
  }
}

// A term of the Hamiltonian took a basis state to a state that is not in the
// basis. With momentum sectors this also happens for states whose orbit has
// zero norm in the sector; those elements are dropped. If instead the term
// changes a quantum number that the basis fixes (CONSERVED_QUANTUMNUMBERS or a
// <QN>_total parameter), dropping them would change the spectrum: stop.
template <class M, class G>
void hamiltonian_matrix<M,G>::check_sector(const state_type& state, const state_type& newstate,
                                           const std::string& term, const std::vector<std::size_t>& sites) const
{
  const typename basis_descriptor_type::constraints_type& constraints = basis_.get_basis().constraints();
  for (std::size_t c=0;c<constraints.size();++c) {
    typename basis_descriptor_type::total_type change;
    for (std::size_t r=0;r<sites.size();++r) {
      const std::size_t s=sites[r];
      const std::size_t k=get_quantumnumber_index(constraints[c].first,basis_.get_site_basis(s));
      if (k==basis_.get_site_basis(s).size())
        continue;
      change += get_quantumnumber(basis_[s][newstate[s]],k) - get_quantumnumber(basis_[s][state[s]],k);
    }
    if (change!=typename basis_descriptor_type::total_type(0)) {
      const std::string& qn=constraints[c].first;
      std::string value;   // the <CONSTRAINT value=...> of the basis, e.g. Nup_total
      const typename basis_descriptor_type::unevaluated_constraints_type& all = basis_.get_basis().all_constraints();
      for (std::size_t a=0;a<all.size();++a)
        if (all[a].first==qn)
          value=boost::lexical_cast<std::string>(all[a].second);
      std::ostringstream msg;
      msg << term << " on site" << (sites.size()>1 ? "s" : "");
      for (std::size_t r=0;r<sites.size();++r)
        msg << " " << sites[r];
      msg << " changes " << qn << ", which the basis fixes in this sector (" << qn << " = "
          << value << " = " << constraints[c].second << "). Remove " << qn
          << " from CONSERVED_QUANTUMNUMBERS and leave " << value << " undefined.";
      boost::throw_exception(std::runtime_error(msg.str()));
    }
  }
}

// Fill one n-site term on n-site group p. The local table (see nsiteoperator.h)
// stores each product reordered as O_0 O_1 ... O_{n-1} (role order). Acting on
// a basis state, O_{n-1} is applied first; each odd O_r picks up the
// Jordan-Wigner sign (-1)^(number of fermions on sites before p.sites[r]),
// counted in the state as modified by the operators already applied. This is
// the rule the bond fill below encodes in closed form for n = 2.
template <class M, class G> template <class STATES, class V, class W>
void hamiltonian_matrix<M,G>::apply_operator(const STATES& states, const nsite_matrix<value_type>& mat, const nsite_group& p, const V& x, W& y) const
{
  if (!built_basis_)                     // the parity table below reads basis_
    build_basis();
  const std::size_t n = p.sites.size();
  const std::size_t nsites = graph_.num_sites();
  // fermion parity of every local state of every site
  std::vector<std::vector<bool> > parity(nsites);
  for (std::size_t s=0;s<nsites;++s)
    for (std::size_t k=0;k<basis_[s].size();++k)
      parity[s].push_back(is_fermionic(model_.model().basis().site_basis(graph_.site_type(s)),basis_[s][k]));

  std::vector<std::size_t> in(n), out(n);
  std::vector<bool> prefix(nsites+1);
  for (std::size_t i=0;i<dimension();++i) {         // loop over source states
    state_type state=states[i];
    for (std::size_t r=0;r<n;++r)
      in[r]=state[p.sites[r]];
    const typename nsite_matrix<value_type>::row_type& row = mat.row(mat.index(in));
    if (row.empty())
      continue;
    prefix[0]=false;                                // parity of sites < s
    for (std::size_t s=0;s<nsites;++s)
      prefix[s+1] = prefix[s] != parity[s][state[s]];
    for (std::size_t e=0;e<row.size();++e) {
      mat.decode(row[e].first,out);
      state_type newstate=state;
      for (std::size_t r=0;r<n;++r)
        newstate[p.sites[r]]=out[r];
      std::complex<double> phase;
      std::size_t j;
      boost::tie(j,phase) = states.index_and_phase(newstate);
      if (j>=dimension()) {
        if (building_hamiltonian_)
          check_sector(state,newstate,"<NSITETERM> of type "+boost::lexical_cast<std::string>(p.type),p.sites);
        continue;
      }
      bool f=false;
      for (std::size_t r=n;r-->0;) {
        std::size_t s=p.sites[r];
        if (parity[s][in[r]]!=parity[s][out[r]]) {  // O_r is odd
          bool before=prefix[s];
          for (std::size_t q=r+1;q<n;++q)           // already applied O_q
            if (p.sites[q]<s && parity[p.sites[q]][in[q]]!=parity[p.sites[q]][out[q]])
              before=!before;
          if (before)
            f=!f;
        }
      }
      value_type val=row[e].second;
      if (f)
        val=-val;
      val *= numeric_cast<value_type>(phase * states.normalization(j)/states.normalization(i));
      y[i] += val*x[j];                   // set matrix element
      simplify(y[i]);               // optionally simplify a symbolic expression
    }
  }
}

template <class M, class G> template <class STATES, class V, class W>
void hamiltonian_matrix<M,G>::apply_operator(const STATES& states, const SiteOperator& op, site_descriptor s, const V& x, W& y) const
{
  multi_array<value_type,2> mat = local_matrix(op,s);
  for (unsigned int i=0;i<dimension();++i) {           // loop basis states
    state_type state=states[i];               // get source state
    int is=state[s];                          // get site basis index
    for (unsigned int js=0;js<basis_[s].size();++js) { // loop over target site states
      value_type val=mat[is][js];                      // get matrix element
      if (numeric::is_nonzero(val)) {   // if matrix element is nonzero
        state_type newstate=state;            // prepare target state
        newstate[s]=js;                       // build target state
        std::complex<double> phase;
        std::size_t j;       
        boost::tie(j,phase) = states.index_and_phase(newstate);       // lookup target state
        if (j<dimension()) {
          val *= numeric_cast<value_type>(phase * states.normalization(j)/states.normalization(i));
          y[i] += val*x[j];                   // set matrix element
          simplify(y[i]);               // optionally simplify a symbolic expression
       }
        else if (building_hamiltonian_)
          check_sector(state,newstate,"<SITETERM>",std::vector<std::size_t>(1,s));
      }
    }
  }
}

template <class M, class G> template <class V, class W>
void hamiltonian_matrix<M,G>::apply_operator(const BondOperator& op, bond_descriptor b, const V& x, W& y) const
{
  apply_operator(local_matrix(op,b),graph_.source(b),graph_.target(b),x,y);
}  


template <class M, class G> template <class V, class W>
void hamiltonian_matrix<M,G>::apply_operator(const BondOperator& op, site_descriptor s1, site_descriptor s2, const V& x, W& y) const
{
  apply_operator(local_matrix(op,s1,s2),s1,s2,x,y);
}


template <class M, class G> template <class STATES, class V, class W>
void hamiltonian_matrix<M,G>::apply_operator(const STATES& states, const boost::multi_array<std::pair<value_type,bool>,4>& mat, site_descriptor s1, site_descriptor s2, const V& x, W& y) const
{
  for (std::size_t i=0;i<dimension();++i) {             // loop over source states
    state_type state=states[i];                   // get source state
    int is1=state[s1];                            // get source site states
    int is2=state[s2];
    for (unsigned int js1=0;js1<basis_[s1].size();++js1) { // loop over target site states
      for (unsigned int js2=0;js2<basis_[s2].size();++js2) {
        value_type val=mat[is1][is2][js1][js2].first;      // get matrix element
        if (numeric::is_nonzero(val)) {              // if nonzero matrix element
          state_type newstate=state;              // prepare target state
          newstate[s1]=js1;                       // build target state
          newstate[s2]=js2;
          std::complex<double> phase;
          std::size_t j;       
          boost::tie(j,phase) = states.index_and_phase(newstate);       // lookup target state
          if (j<dimension()) {
            if (mat[is1][is2][js1][js2].second) {
              // calculate fermionic sign
              bool f=(s2>=s1);
              int start = (std::min)(s1,s2);
              int end = (std::max)(s1,s2);

              for (int i=start;i<end;++i)
                if (is_fermionic(model_.model().basis().site_basis(graph_.site_type(i)),basis_[i][state[i]]))
                  f=!f;
              if (f)
                val=-val;
            }
            val *= numeric_cast<value_type>(phase * states.normalization(j)/states.normalization(i));
            y[i] += val*x[j];                   // set matrix element
            simplify(y[i]);               // optionally simplify a symbolic expression
          }
          else if (building_hamiltonian_) {
            std::vector<std::size_t> sites;
            sites.push_back(s1);
            sites.push_back(s2);
            check_sector(state,newstate,"<BONDTERM>",sites);
          }
        }
      }
    }
  }
}



template <class M, class G> template <class V, class W>
void hamiltonian_matrix<M,G>::apply_operator(const std::string& op, site_descriptor s, const V& x, W& y) const
{
  if (model_.has_site_operator(op))
    apply_operator(model_.get_site_operator(op),s,x,y);
  else
    apply_operator(SiteOperator(op+"(i)","i"),s,x,y);
}


template <class M, class G> template <class V, class W>
void hamiltonian_matrix<M,G>::apply_operator(const std::string& name, const V& x, W& y) const
{ 
  if (model_.has_site_operator(name))
    apply_operator(model_.get_site_operator(name),x,y); 
  else if (model_.has_bond_operator(name))
    apply_operator(model_.get_bond_operator(name),x,y);
  else if (model_.has_global_operator(name))
    apply_operator(model_.get_global_operator(name),x,y);
  else // assume site operator
    apply_operator(SiteOperator(name+"(i)","i"),x,y);
}


template <class M, class G> template <class V, class W>
void hamiltonian_matrix<M,G>::apply_operator(const SiteOperator& op, const V& x, W& y) const
{
  for (site_iterator it=graph_.sites().first; it!=graph_.sites().second ; ++it) 
    apply_operator(op,*it,x,y);
}
  
template <class M, class G> template <class V, class W>
void hamiltonian_matrix<M,G>::apply_operator(const BondOperator& op, const V& x, W& y) const
{
  for (bond_iterator it=graph_.bonds().first; it!=graph_.bonds().second ; ++it) 
    apply_operator(op,*it,x,y);
}

template <class M, class G> template <class V, class W>
void hamiltonian_matrix<M,G>::apply_operator(const GlobalOperator& op, const V& x, W& y) const
{
  // loop over sites
  for (site_iterator it=graph_.sites().first; it!=graph_.sites().second ; ++it)
    apply_operator(op.site_term(graph_.site_type(*it)),*it,x,y);

  // loop over bonds
  for (bond_iterator it=graph_.bonds().first; it!=graph_.bonds().second ; ++it) 
      apply_operator(op.bond_term(graph_.bond_type(*it)),*it,x,y);

  // loop over n-site groups (n-site terms). The local table depends only on the
  // n-site group type and the site types, so build it once per combination: the
  // expression evaluator is slow.
  if (graph_.nsite_groups().empty() || !op.has_nsite_terms())
    return;
  typedef std::pair<unsigned int,std::vector<unsigned int> > key_type;
  std::map<key_type,nsite_matrix<value_type> > cache;
  for (std::size_t q=0;q<graph_.nsite_groups().size();++q) {
    const nsite_group& p = graph_.nsite_groups()[q];
    key_type key(p.type,std::vector<unsigned int>());
    for (std::size_t r=0;r<p.sites.size();++r)
      key.second.push_back(graph_.site_type(p.sites[r]));
    typename std::map<key_type,nsite_matrix<value_type> >::iterator it = cache.find(key);
    if (it==cache.end()) {
      NSiteOperator term = op.nsite_term(p.type);
      nsite_matrix<value_type> mat;
      if (!term.term().empty()) {
        mat = local_matrix(term,p);
        warn_if_cutoff(mat,p);
      }
      it = cache.insert(std::make_pair(key,mat)).first;
    }
    if (it->second.size())
      apply_operator(it->second,p,x,y);
  }
}




template <class M, class G>
void hamiltonian_matrix<M,G>::build_basis() const
{
  basis_descriptor_type b = model_.basis();
  b.set_parameters(parms);
  basis_ = basis_states_descriptor<short>(b,graph_.graph());
  if (uses_translation_invariance()) {
    std::vector<Expression> k;
    read_vector_resize(parms["TOTAL_MOMENTUM"],k);
    ParameterEvaluator eval(parms);
    vector_type total_momentum;
    for (unsigned i=0;i<k.size();++i)
      total_momentum.push_back(std::real(k[i].value(eval)));
    bloch_states = bloch_basis_states_type(basis_,graph_.translations(total_momentum));
  }
  else
    states = basis_states_type(basis_);
  built_basis_ = true;
}    

template <class M, class G>
void hamiltonian_matrix<M,G>::build() const
{
  // build matrix
  if (!built_basis_)
    build_basis();
  Disorder::seed(parms.value_or_default("DISORDER_SEED",0));
  matrix_ = detail::initialized_matrix(matrix_type(),dimension(),dimension());
  built_matrix_ = true;
  building_hamiltonian_ = true;
  try {
    add_operator_matrix(matrix_,model_.model());
  }
  catch (...) {
    building_hamiltonian_ = false;
    built_matrix_ = false;
    throw;
  }
  building_hamiltonian_ = false;

  //std::cerr << "Time to build matrix: " << t.elapsed() << "\n";
/*
  if (this->parms.value_or_default("CHECK_SYMMETRY",false)) {
    std::cerr << "Checking symmetry\n";
    for (unsigned i=0;i<dimension();++i)
      for (unsigned j=0;j<i;++j)
        if (std::abs(static_cast<value_type>(matrix_(i,j))-numeric::conj(static_cast<value_type>(matrix_(j,i)))) > 1e-10) {
          std::cerr << "Symmetry problem: " << i << " " << j << " " << static_cast<value_type>(matrix_(i,j)) << " " << numeric::conj(static_cast<value_type>(matrix_(j,i))) << "\n";
          std::cout << basis() << "\n";
          std::cout << states << "\n";
          std::cout << bloch_states << "\n";
          std::cout << matrix() << "\n";
          std::abort();
        }
    std::cerr << "Checked symmetry\n";
  }
*/
}

} // end name space alps

#endif // ALPS_MODEL_MATRIX_HPP

