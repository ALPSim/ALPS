/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2005 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

/* $Id$ */

#ifndef ALPS_EXPRESSION_EVALUATABLE_H
#define ALPS_EXPRESSION_EVALUATABLE_H

#include <alps/expression/expression_fwd.h>

namespace alps {
namespace expression {


template<class T>
class Evaluatable {
public:
  typedef T value_type;

  Evaluatable() {}
  virtual ~Evaluatable() {}
  virtual value_type value(const Evaluator<T>& =Evaluator<T>(), bool=false) const = 0;
  virtual bool can_evaluate(const Evaluator<T>& =Evaluator<T>(), bool=false) const = 0;
  virtual void output(std::ostream&) const = 0;
  virtual Evaluatable* clone() const = 0;
  virtual boost::shared_ptr<Evaluatable> flatten_one() { return boost::shared_ptr<Evaluatable>(); }
  virtual Evaluatable* partial_evaluate_replace(const Evaluator<T>& =Evaluator<T>(),bool=false) { return this; }
  virtual bool is_single_term() const { return false; }
  virtual Term<T> term() const;
  virtual bool depends_on(const std::string&) const { return false; }
};


template<class T>
inline Term<T> Evaluatable<T>::term() const { return Term<T>(); }

} // end namespace expression
} // end namespace alps

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
namespace alps {
namespace expression {
#endif

template<class T>
inline std::ostream& operator<<(std::ostream& os, const alps::expression::Evaluatable<T>& e)
{
  e.output(os);
  return os;
}

#ifndef BOOST_NO_OPERATORS_IN_NAMESPACE
} // end namespace expression
} // end namespace alps
#endif

#endif // ! ALPS_EXPRESSION_IMPL_H
