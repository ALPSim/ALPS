/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2005 by Matthias Troyer <troyer@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the “Software”),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
* OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
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
