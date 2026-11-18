/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2006 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_MODEL_OPERATOR_H
#define ALPS_MODEL_OPERATOR_H

#include <alps/expression.h>
#include <alps/parameter.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace alps {

// An expression, term or factor as it prints.
template <class X>
std::string expression_text(const X& x)
{
  std::ostringstream out;
  out << x;
  return out.str();
}

// Throws if a site operator on one of the given sites stands in a denominator,
// e.g. 1/Sz(i). The operator evaluators apply such an operator as if it were in
// the numerator, and simplify() drops the division (1/Sz(i) becomes Sz(i)), so
// call this before the term is evaluated or simplified.
template <class T>
void refuse_operator_in_denominator(const expression::Expression<T>& term, const std::vector<std::string>& sites)
{
  for (typename expression::Expression<T>::term_iterator t = term.terms().first; t != term.terms().second; ++t)
    for (typename expression::Term<T>::factor_iterator f = t->factors().first; f != t->factors().second; ++f) {
      bool on_site = false;
      for (std::size_t s = 0; s < sites.size(); ++s)
        on_site = on_site || f->depends_on(sites[s]);
      if (!on_site)
        continue;
      if (f->is_inverse())
        throw std::runtime_error("Cannot divide by a site operator: "
          + expression_text(*t) + " divides by " + expression_text(*f));
      // look inside parentheses, also under a power: K*(Sz(j)/Sz(i)), (1/Sz(i))^2
      if (const expression::Expression<T>* inner = dynamic_cast<const expression::Expression<T>*>(f->evaluatable()))
        refuse_operator_in_denominator(*inner, sites);
    }
}

template <class T=std::complex<double> >
class OperatorEvaluator : public expression::ParameterEvaluator<T>
{
public:
  typedef expression::ParameterEvaluator<T> super_type;
  typedef typename super_type::value_type value_type;
  
  OperatorEvaluator(const Parameters& p)
    : super_type(p) {}
  typename super_type::Direction direction() const { return super_type::right_to_left; }

  value_type evaluate(const std::string& name, bool isarg=false) const
  { return super_type::partial_evaluate(name,isarg).value();}

  value_type evaluate_function(const std::string& name, const expression::Expression<T>& arg,bool isarg=false) const
  { return super_type::partial_evaluate_function(name,arg,isarg).value();}

  value_type evaluate_function(const std::string& name, const std::vector<expression::Expression<T> >& args,bool isarg=false) const
  { return super_type::partial_evaluate_function(name,args,isarg).value();}
};

} // namespace alps

#endif
