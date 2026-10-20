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

// Terms acting on n sites (n arbitrary), e.g. a 4-site term
//
//   <NSITETERM type="0" sites="i j k l">
//     K*(cdag_up(l)*cdag_down(k)*c_down(j)*c_up(i) + h.c.)
//   </NSITETERM>
//
// applied to every <NSITE> of the lattice. Modeled on BondOperator /
// BondOperatorSplitter (bondoperator.h, bondterm.h), with the pair of sites
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
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace alps {

class ModelLibrary;

// Local table of an n-site term: for each input tuple of local states the
// nonzero (output tuple, value) pairs. Tuples are mixed-radix indices, role 0
// most significant.
//
// The table is stored factorized: the term is a sum of products
// coef * O_0 * O_1 * ... * O_{n-1}, one single-site operator per role (or the
// identity). A row is assembled from these the first time it is asked for, so
// only the rows that basis states actually reach are built; the full table
// has prod(dims) rows, which for large n-site groups or boson cutoffs is far more
// than any one sector uses. Not thread safe.
template <class T>
struct nsite_matrix
{
  typedef std::vector<std::pair<std::size_t, T> > row_type;
  // a single-site operator: for each input local state, the (output, value) pairs
  typedef std::vector<row_type> local_type;
  struct product_type
  {
    T coef;
    std::vector<int> ops;   // per role: index into locals, or -1 for the identity
  };

  std::vector<std::size_t> dims;
  std::vector<local_type> locals;
  std::vector<product_type> products;
  // quantum numbers whose sum over the sites some product changes
  std::set<std::string> nonconserved;

  // number of rows (input tuples); 0 if the term has no nonzero product
  std::size_t size() const
  {
    if (products.empty())
      return 0;
    std::size_t total = 1;
    for (std::size_t r = 0; r < dims.size(); ++r)
      total *= dims[r];
    return total;
  }

  const row_type& row(std::size_t idx) const
  {
    typename std::unordered_map<std::size_t, row_type>::const_iterator it = rows_.find(idx);
    if (it != rows_.end())
      return it->second;
    std::vector<std::size_t> in(dims.size()), out(dims.size());
    decode(idx, in);
    std::map<std::size_t, T> acc;
    for (std::size_t p = 0; p < products.size(); ++p)
      expand(products[p], 0, products[p].coef, in, out, acc);
    row_type& result = rows_[idx];
    for (typename std::map<std::size_t, T>::iterator e = acc.begin(); e != acc.end(); ++e) {
      simplify(e->second);
      if (numeric::is_nonzero(e->second))
        result.push_back(*e);
    }
    return result;
  }

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

private:
  void expand(const product_type& p, std::size_t r, const T& value, const std::vector<std::size_t>& in,
              std::vector<std::size_t>& out, std::map<std::size_t, T>& acc) const
  {
    if (r == dims.size()) {
      acc[index(out)] += value;
      return;
    }
    if (p.ops[r] < 0) {
      out[r] = in[r];
      expand(p, r + 1, value, in, out, acc);
      return;
    }
    const row_type& column = locals[p.ops[r]][in[r]];
    for (std::size_t e = 0; e < column.size(); ++e) {
      out[r] = column[e].first;
      expand(p, r + 1, value * column[e].second, in, out, acc);
    }
  }

  mutable std::unordered_map<std::size_t, row_type> rows_;
};


// Splits one product term of an n-site operator into a coefficient and one
// string of site operators per role, like BondOperatorSplitter for two sites.
template <class I, class T = std::complex<double> >
class NSiteOperatorSplitter : public OperatorEvaluator<T>
{
private:
  typedef OperatorEvaluator<T> super_type;

public:
  NSiteOperatorSplitter(const std::vector<SiteBasisDescriptor<I> >& bases,
                        const std::vector<std::string>& sites,
                        const Parameters& p)
    : super_type(p), bases_(bases), sites_(sites), ops_(sites.size()),
      used_(sites.size(), false), odd_(sites.size(), false) {}

  bool can_evaluate_function(const std::string& name, const expression::Expression<T>& arg, bool isarg = false) const
  {
    return role(arg) >= 0 || expression::ParameterEvaluator<T>::can_evaluate_function(name, arg, isarg);
  }
  expression::Expression<T> partial_evaluate_function(const std::string& name,
                                                      const expression::Expression<T>& argument, bool = false) const;

  // the operators of the term on role r, in their original order
  const expression::Term<T>& site_operators(std::size_t r) const { return ops_[r]; }
  bool used(std::size_t r) const { return used_[r]; }
  // the term has an odd number of fermionic operators
  bool odd() const
  {
    bool f = false;
    for (std::size_t r = 0; r < odd_.size(); ++r)
      if (odd_[r])
        f = !f;
    return f;
  }

private:
  int role(const expression::Expression<T>& arg) const
  {
    for (std::size_t r = 0; r < sites_.size(); ++r)
      if (arg == sites_[r])
        return r;
    return -1;
  }

  const std::vector<SiteBasisDescriptor<I> >& bases_;
  const std::vector<std::string>& sites_;
  mutable std::vector<expression::Term<T> > ops_;
  mutable std::vector<bool> used_;
  mutable std::vector<bool> odd_;   // parity of the fermionic operators seen on each role
};

// Operators are applied right to left (OperatorEvaluator::direction), so all
// operators already seen stand to the right of the current one. Moving the
// current operator on role r into role order means moving it right past the
// already-seen operators on roles q < r: one sign per odd one among them.
// For n = 2 this is BondOperatorSplitter's rule with the roles swapped.
template <class I, class T>
expression::Expression<T> NSiteOperatorSplitter<I, T>::partial_evaluate_function(const std::string& name, const expression::Expression<T>& arg, bool isarg) const
{
  int r = role(arg);
  if (r >= 0) {
    // Replaced by its sign, an operator inside another function such as
    // sqrt(n(i)) would enter the table as n(i) itself.
    if (isarg)
      boost::throw_exception(std::runtime_error("n-site term: the site operator " + name + "("
        + boost::lexical_cast<std::string>(arg) + ") is the argument of another function, which is not supported"));
    ops_[r] = expression::Function<T>(name, arg) * ops_[r];
    used_[r] = true;
    double sign = 1.;
    if (bases_[r].is_fermionic(name)) {
      for (int q = 0; q < r; ++q)
        if (odd_[q])
          sign = -sign;
      odd_[r] = !odd_[r];
    }
    return expression::Expression<T>(sign);
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
  for (std::size_t r = 0; r < n; ++r) {
    states.push_back(site_basis<I>(bases[r]));
    mat.dims.push_back(states[r].size());
  }
  // rows are keyed by the mixed-radix index of the local states (and size()
  // is its range), so prod(dims) must fit in a std::size_t
  std::size_t configurations = 1;
  for (std::size_t r = 0; r < n; ++r) {
    if (mat.dims[r] && configurations > std::numeric_limits<std::size_t>::max() / mat.dims[r]) {
      std::string d;
      for (std::size_t q = 0; q < n; ++q)
        d += (q ? " x " : "") + boost::lexical_cast<std::string>(mat.dims[q]);
      boost::throw_exception(std::runtime_error("n-site term on " + boost::lexical_cast<std::string>(n)
        + " sites: the number of local configurations (" + d + ") does not fit in a "
        + boost::lexical_cast<std::string>(std::numeric_limits<std::size_t>::digits)
        + "-bit index. Write the term on n-site groups with fewer sites."));
    }
    configurations *= mat.dims[r];
  }

  // parse expression and store it as sum of terms
  expression::Expression<value_type> ex(term());
  ex.flatten();
  ex.simplify();

  // single-site operators already built, by (role, operator string)
  std::map<std::pair<std::size_t, std::string>, int> built;
  for (typename expression::Expression<value_type>::term_iterator tit = ex.terms().first; tit != ex.terms().second; ++tit) {
    NSiteOperatorSplitter<I, value_type> splitter(bases, sites(), parms);
    expression::Term<value_type> term(*tit);
    // The splitter replaces each site operator by a sign, so under a power
    // the operator would be applied once and only the sign raised.
    for (typename expression::Term<value_type>::factor_iterator f = term.factors().first; f != term.factors().second; ++f)
      if (!f->unit_power())
        for (std::size_t r = 0; r < n; ++r)
          if (f->depends_on(sites()[r]))
            boost::throw_exception(std::runtime_error("n-site term: powers of site operators are not supported in "
              + boost::lexical_cast<std::string>(*f) + "; write the product out, e.g. n(i)*n(i) for n(i)^2"));
    term.partial_evaluate(splitter);
    term.simplify();
    if (!numeric::is_nonzero(term))
      continue;
    if (boost::is_arithmetic<T>::value || is_complex<T>::value)
      if (!can_evaluate(boost::lexical_cast<std::string>(term)))
        boost::throw_exception(std::runtime_error("Cannot evaluate expression " + boost::lexical_cast<std::string>(term)));
    typename nsite_matrix<T>::product_type prod;
    prod.coef = evaluate<T>(term);
    if (!numeric::is_nonzero(prod.coef))
      continue;

    bool zero = false;
    std::map<std::string, half_integer<int> > change;   // total over the roles
    for (std::size_t r = 0; r < n && !zero; ++r) {
      if (!splitter.used(r)) {
        prod.ops.push_back(-1);
        continue;
      }
      std::pair<std::size_t, std::string> key(r, boost::lexical_cast<std::string>(splitter.site_operators(r)));
      std::map<std::pair<std::size_t, std::string>, int>::const_iterator it = built.find(key);
      if (it == built.end()) {
        multi_array<std::pair<T, bool>, 2> m = SiteOperator(key.second, sites()[r]).template matrix<T, I>(bases[r], parms);
        typename nsite_matrix<T>::local_type local(mat.dims[r]);
        for (std::size_t i = 0; i < mat.dims[r]; ++i)
          for (std::size_t j = 0; j < mat.dims[r]; ++j)
            if (numeric::is_nonzero(m[i][j].first))
              local[i].push_back(std::make_pair(j, m[i][j].first));
        mat.locals.push_back(local);
        it = built.insert(std::make_pair(key, int(mat.locals.size()) - 1)).first;
      }
      prod.ops.push_back(it->second);
      // the change of each quantum number on this role, read off one nonzero element
      const typename nsite_matrix<T>::local_type& local = mat.locals[it->second];
      zero = true;
      for (std::size_t i = 0; i < local.size() && zero; ++i)
        if (!local[i].empty()) {
          zero = false;
          for (std::size_t q = 0; q < bases[r].size(); ++q)
            change[bases[r][q].name()] += get_quantumnumber(states[r][local[i][0].first], q)
                                          - get_quantumnumber(states[r][i], q);
        }
    }
    if (zero)
      continue;
    if (splitter.odd())
      boost::throw_exception(std::runtime_error("n-site term contains an odd number of fermion operators: "
                                                + boost::lexical_cast<std::string>(*tit)));
    for (std::map<std::string, half_integer<int> >::const_iterator c = change.begin(); c != change.end(); ++c)
      if (c->second != half_integer<int>(0))
        mat.nonconserved.insert(c->first);
    mat.products.push_back(prod);
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
