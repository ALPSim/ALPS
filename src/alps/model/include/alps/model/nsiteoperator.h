/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#ifndef ALPS_MODEL_NSITEOPERATOR_H
#define ALPS_MODEL_NSITEOPERATOR_H

// Terms acting on n sites (n arbitrary), e.g. a 4-site plaquette term
//
//   <PLAQUETTETERM type="0" sites="i j k l">
//     K*(cdag_up(l)*cdag_down(k)*c_down(j)*c_up(i) + h.c.)
//   </PLAQUETTETERM>
//
// applied to every <PLAQUETTE> of the lattice. Modeled on BondOperator /
// BondOperatorEvaluator (bondoperator.h, bondterm.h), with the pair of sites
// replaced by an ordered list of roles.
//
// Fermion convention of the local table. Each product of operators in the term
// is reordered so that all operators on role 0 stand leftmost, then role 1, ...
// (operators on the same role keep their order): O_0 O_1 ... O_{n-1}. The sign
// of that reordering is folded into the table. Applying the table to a basis
// state then needs only the usual Jordan-Wigner string of each odd O_r, which
// hamiltonian_matrix computes. (BondOperator uses the opposite role order,
// O_target O_source; both are self-consistent.)

#include <alps/config.h>
#include <alps/model/operator.h>
#include <alps/model/siteoperator.h>
#include <alps/model/siteterm.h>
#include <alps/model/sitestate.h>
#include <alps/model/sitebasisstates.h>
#include <alps/numeric/is_nonzero.hpp>
#include <alps/type_traits/is_complex.hpp>
#include <alps/expression.h>
#include <alps/parameter.h>
#include <boost/lexical_cast.hpp>
#include <boost/throw_exception.hpp>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace alps {

class ModelLibrary;

// Sparse local matrix of an n-site term: for each input tuple of local states
// the nonzero (output tuple, value) pairs. Tuples are mixed-radix indices,
// role 0 most significant. A dense [dim]^(2n) array would be mostly zeros
// (n=4, dim=4: 65536 entries).
template <class T>
struct nsite_matrix
{
  typedef std::vector<std::pair<std::size_t, T> > row_type;
  std::vector<std::size_t> dims;
  std::vector<row_type> rows;
  // quantum numbers whose sum over the sites some nonzero entry changes
  std::set<std::string> nonconserved;

  std::size_t size() const { return rows.size(); }

  template <class V>
  std::size_t index(const V& local) const
  {
    std::size_t idx = 0;
    for (std::size_t r = 0; r < dims.size(); ++r)
      idx = idx * dims[r] + local[r];
    return idx;
  }

  template <class V>
  void decode(std::size_t idx, V& local) const
  {
    for (std::size_t r = dims.size(); r-- > 0;) {
      local[r] = idx % dims[r];
      idx /= dims[r];
    }
  }
};


template <class I, class T = std::complex<double>, class STATE = site_state<I> >
class NSiteOperatorEvaluator : public OperatorEvaluator<T>
{
private:
  typedef OperatorEvaluator<T> super_type;
  typedef NSiteOperatorEvaluator<I, T, STATE> SELF_;

public:
  NSiteOperatorEvaluator(const std::vector<STATE>& states,
                         const std::vector<SiteBasisDescriptor<I> >& bases,
                         const std::vector<std::string>& sites,
                         const Parameters& p)
    : super_type(p)
  {
    for (std::size_t r = 0; r < sites.size(); ++r)
      site_.push_back(SiteOperatorEvaluator<I, T, STATE>(states[r], bases[r], p, sites[r]));
  }

  bool can_evaluate_function(const std::string&, const expression::Expression<T>&, bool = false) const;
  expression::Expression<T> partial_evaluate_function(const std::string& name,
                                                      const expression::Expression<T>& argument, bool = false) const;

  std::size_t size() const { return site_.size(); }
  const STATE& state(std::size_t r) const { return site_[r].state(); }
  // parity of the number of fermionic operators applied on role r
  bool fermionic(std::size_t r) const { return site_[r].fermionic(); }

  bool has_operator(const std::string& name, const expression::Expression<T>& arg) const
  {
    int r = role(arg);
    return r >= 0 && site_[r].has_operator(name);
  }

private:
  int role(const expression::Expression<T>& arg) const
  {
    for (std::size_t r = 0; r < site_.size(); ++r)
      if (arg == site_[r].site())
        return r;
    return -1;
  }

  std::vector<SiteOperatorEvaluator<I, T, STATE> > site_;
};


template <class I, class T, class STATE>
bool NSiteOperatorEvaluator<I, T, STATE>::can_evaluate_function(const std::string& name, const expression::Expression<T>& arg, bool isarg) const
{
  if (has_operator(name, arg)) {
    SELF_ eval(*this);
    return eval.partial_evaluate_function(name, arg, isarg).can_evaluate(expression::ParameterEvaluator<T>(*this), isarg);
  }
  else
    return expression::ParameterEvaluator<T>::can_evaluate_function(name, arg, isarg);
}

// Operators are applied right to left (OperatorEvaluator::direction), so all
// operators already applied stand to the right of the current one. Moving the
// current operator on role r into role order means moving it right past the
// already-applied operators on roles q < r: one sign per odd one among them.
// For n = 2 this is BondOperatorEvaluator's rule with the roles swapped.
template <class I, class T, class STATE>
expression::Expression<T> NSiteOperatorEvaluator<I, T, STATE>::partial_evaluate_function(const std::string& name, const expression::Expression<T>& arg, bool isarg) const
{
  if (has_operator(name, arg)) {
    int r = role(arg);
    bool f = site_[r].fermionic();
    expression::Expression<T> e = site_[r].partial_evaluate_function(name, site_[r].site(), isarg);
    if (f != site_[r].fermionic() && numeric::is_nonzero(e)) {
      bool odd = false;
      for (int q = 0; q < r; ++q)
        if (site_[q].fermionic())
          odd = !odd;
      if (odd)
        e.negate();
    }
    return e;
  }
  return expression::ParameterEvaluator<T>(*this).partial_evaluate_function(name, arg, isarg);
}


class ALPS_DECL NSiteOperator
{
public:
  NSiteOperator() {}
  NSiteOperator(const std::string& term, const std::vector<std::string>& sites)
    : term_(term), sites_(sites) {}
  NSiteOperator(const XMLTag& tag, std::istream& in) { read_xml(tag, in); }
  NSiteOperator(NSiteOperator const& op, std::string const& t, Parameters const& p)
    : name_(op.name_), term_(t), sites_(op.sites_), parms_(p) {}

  void read_xml(const XMLTag&, std::istream&);
  void write_xml(oxstream&) const;

  const std::string& name() const { return name_; }
  const std::string& term() const { return term_; }
  std::string& term() { return term_; }
  const std::vector<std::string>& sites() const { return sites_; }
  std::size_t num_sites() const { return sites_.size(); }
  Parameters const& parms() const { return parms_; }

  void substitute_operators(const ModelLibrary& m, const Parameters& p = Parameters());

  // Local table on sites with the given site bases (one per role).
  template <class T, class I>
  nsite_matrix<T> matrix(const std::vector<SiteBasisDescriptor<I> >& bases,
                         const Parameters& p = Parameters()) const;

protected:
  std::string sites_string() const;

private:
  std::string name_;
  std::string term_;
  std::vector<std::string> sites_;
  Parameters parms_;
};


class ALPS_DECL NSiteTermDescriptor : public NSiteOperator
{
public:
  NSiteTermDescriptor() : NSiteOperator(), type_(-2) {}
  NSiteTermDescriptor(const std::string& term, const std::vector<std::string>& sites)
    : NSiteOperator(term, sites), type_(-2) {}
  NSiteTermDescriptor(const XMLTag&, std::istream&);
  NSiteTermDescriptor(NSiteTermDescriptor const& t, std::string const& term, Parameters const& p, unsigned int type)
    : NSiteOperator(t, term, p), type_(type) {}

  const NSiteOperator& nsite_operator() const { return static_cast<const NSiteOperator&>(*this); }
  void write_xml(oxstream&) const;

  bool match_type(int type) const { return type_ == -1 || type == type_; }
  int type() const { return type_; }
private:
  int type_;
};


template <class T, class I>
nsite_matrix<T> NSiteOperator::matrix(const std::vector<SiteBasisDescriptor<I> >& b, const Parameters& p) const
{
  typedef typename expression_value_type_traits<T>::value_type value_type;

  const std::size_t n = num_sites();
  if (b.size() != n)
    boost::throw_exception(std::runtime_error("n-site term with " + boost::lexical_cast<std::string>(n)
      + " sites applied to a group of " + boost::lexical_cast<std::string>(b.size()) + " sites"));

  std::vector<SiteBasisDescriptor<I> > bases(b);
  Parameters parms(p);
  for (std::size_t r = 0; r < n; ++r)
    bases[r].set_parameters(p);
  for (std::size_t r = 0; r < n; ++r)
    parms.copy_undefined(bases[r].get_parameters());

  std::vector<site_basis<I> > states;
  nsite_matrix<T> mat;
  std::size_t total = 1;
  for (std::size_t r = 0; r < n; ++r) {
    states.push_back(site_basis<I>(bases[r]));
    mat.dims.push_back(states[r].size());
    total *= states[r].size();
  }
  mat.rows.resize(total);

  // parse expression and store it as sum of terms
  expression::Expression<value_type> ex(term());
  ex.flatten();
  ex.simplify();

  std::vector<std::size_t> in(n);
  std::vector<std::size_t> out(n);
  std::vector<site_state<I> > instates(n);
  for (std::size_t i = 0; i < total; ++i) {
    mat.decode(i, in);
    for (std::size_t r = 0; r < n; ++r)
      instates[r] = states[r][in[r]];
    std::map<std::size_t, T> row;
    for (typename expression::Expression<value_type>::term_iterator tit = ex.terms().first; tit != ex.terms().second; ++tit) {
      NSiteOperatorEvaluator<I, value_type> evaluator(instates, bases, sites(), parms);
      expression::Term<value_type> term(*tit);
      term.partial_evaluate(evaluator);
      if (!numeric::is_nonzero(term))
        continue;
      bool valid = true;
      for (std::size_t r = 0; r < n; ++r) {
        out[r] = states[r].index(evaluator.state(r));
        if (out[r] >= mat.dims[r])
          valid = false;
      }
      if (!valid)
        continue;
      if (boost::is_arithmetic<T>::value || is_complex<T>::value)
        if (!can_evaluate(boost::lexical_cast<std::string>(term)))
          boost::throw_exception(std::runtime_error("Cannot evaluate expression " + boost::lexical_cast<std::string>(term)));
      T value = evaluate<T>(term);
      // An operator that annihilates the state leaves a zero factor which
      // partial_evaluate does not always fold away ("1 * (0)"), and a later
      // operator can make the state valid again (c then cdag on an empty
      // mode). Drop such zeros before looking at the fermion parity.
      if (!numeric::is_nonzero(value))
        continue;
      bool odd = false;
      for (std::size_t r = 0; r < n; ++r)
        if (evaluator.fermionic(r))
          odd = !odd;
      if (odd)
        boost::throw_exception(std::runtime_error("n-site term contains an odd number of fermion operators: "
                                                  + boost::lexical_cast<std::string>(*tit)));
      T& entry = row[mat.index(out)];
      entry += value;
      simplify(entry);
    }
    for (typename std::map<std::size_t, T>::const_iterator it = row.begin(); it != row.end(); ++it)
      if (numeric::is_nonzero(it->second)) {
        mat.rows[i].push_back(*it);
        // record quantum numbers (by name, summed over the roles) it changes
        mat.decode(it->first, out);
        std::map<std::string, half_integer<I> > change;
        for (std::size_t r = 0; r < n; ++r)
          for (std::size_t q = 0; q < bases[r].size(); ++q)
            change[bases[r][q].name()] += get_quantumnumber(states[r][out[r]], q) - get_quantumnumber(instates[r], q);
        for (typename std::map<std::string, half_integer<I> >::const_iterator c = change.begin(); c != change.end(); ++c)
          if (c->second != half_integer<I>(0))
            mat.nonconserved.insert(c->first);
      }
  }
  return mat;
}

} // namespace alps

namespace alps {

inline alps::oxstream& operator<<(alps::oxstream& out, const alps::NSiteTermDescriptor& q)
{
  q.write_xml(out);
  return out;
}

inline std::ostream& operator<<(std::ostream& out, const alps::NSiteTermDescriptor& q)
{
  alps::oxstream xml(out);
  xml << q;
  return out;
}

} // namespace alps

#endif
