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

#ifndef ALPS_EXPRESSION_TRAITS_H
#define ALPS_EXPRESSION_TRAITS_H

#include <alps/expression/expression_fwd.h>

namespace alps {
namespace expression {

template<class T>
struct expression {
  typedef std::complex<T> value_type;
  typedef Expression<value_type> type;
  typedef Term<value_type> term_type;
};

template<class T>
struct expression<std::complex<T> > {
  typedef std::complex<T> value_type;
  typedef Expression<value_type> type;
  typedef Term<value_type> term_type;
};

template<class T>
struct expression<Expression<T> > {
  typedef T value_type;
  typedef Expression<value_type> type;
  typedef Term<value_type> term_type;
};

}


template <class T>
struct expression_value_type_traits {
  typedef typename expression::expression<T>::value_type value_type;
};

//
// function is_zero and is_nonzero
//

template<class T>
bool is_zero(const expression::Expression<T>& x)
{
  std::string s = boost::lexical_cast<std::string>(x);
  return s=="" || s=="0" || s=="0." || s=="-0" || s=="-0.";
}

template<class T>
bool is_zero(const expression::Term<T>& x)
{
  std::string s = boost::lexical_cast<std::string>(x);
  return s=="" || s=="0" || s=="0.";
}


} // end namespace alps

#endif // ! ALPS_EXPRESSION_H
