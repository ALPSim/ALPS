/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2010 by Matthias Troyer <troyer@comp-phys.org>,
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

#ifndef ALPS_EXPRESSION_NUMBER_H
#define ALPS_EXPRESSION_NUMBER_H

#include <alps/expression/expression_fwd.h>
#include <alps/expression/evaluate_helper.h>
#include <alps/type_traits/real_type.hpp>
#include <boost/call_traits.hpp>

namespace alps {
namespace expression {

template<class T>
class Number : public Evaluatable<T> {
public:
  typedef T value_type;
  typedef typename alps::real_type<T>::type real_type;

  Number(typename boost::call_traits<value_type>::param_type x) : val_(x) {}
  value_type value(const Evaluator<T>& =Evaluator<T>(), bool=false) const;
  bool can_evaluate(const Evaluator<T>& =Evaluator<T>(), bool=false) const { return true; }
  void output(std::ostream&) const;
  Evaluatable<T>* clone() const { return new Number<T>(*this); }
private:
  value_type val_;
};

template<class T>
typename Number<T>::value_type Number<T>::value(const Evaluator<T>&, bool) const
{
  return val_;
}

template<class T>
void Number<T>::output(std::ostream& os) const
{
  if (evaluate_helper<T>::imag(val_) == 0)
    os << evaluate_helper<T>::real(val_);
  else
    os << val_;
}

} // end namespace expression
} // end namespace alps

#endif // ! ALPS_EXPRESSION_IMPL_H
