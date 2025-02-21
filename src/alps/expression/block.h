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

#ifndef ALPS_EXPRESSION_BLOCK_H
#define ALPS_EXPRESSION_BLOCK_H

#include <alps/expression/expression_fwd.h>
#include <alps/expression/evaluatable.h>

namespace alps {
namespace expression {

template<class T>
class Block : public Expression<T> {
private:
  typedef Expression<T> BASE_;

public:
  Block(std::istream&);
  Block(const Expression<T>& e) : BASE_(e) {}
  void output(std::ostream&) const;
  Evaluatable<T>* clone() const { return new Block<T>(*this); }
  void flatten();
  boost::shared_ptr<Evaluatable<T> > flatten_one();
  Evaluatable<T>* partial_evaluate_replace(const Evaluator<T>& =Evaluator<T>(),bool=false);
};

//
// implementation of Block<T>
//

template<class T>
Block<T>::Block(std::istream& in) : Expression<T>(in)
{
  char c;
  in >> c;
  if (c != ')' && c != ',')
    boost::throw_exception(std::runtime_error(") or , expected in expression"));
  if (c == ',') {
    // read imaginary part
    Expression<T> ex(in);
    Block<T> bl(ex);
    Term<T> term(bl);
    term *= "I";
    *this += term;
    check_character(in,')',") expected in expression");
  }
}

template<class T>
boost::shared_ptr<Evaluatable<T> > Block<T>::flatten_one()
{
  boost::shared_ptr<Expression<T> > ex = BASE_::flatten_one_expression();
  if (ex)
    return boost::shared_ptr<Evaluatable<T> >(new Block<T>(*ex));
  else
    return boost::shared_ptr<Evaluatable<T> >();
}

template<class T>
void Block<T>::output(std::ostream& os) const
{
  os << "(";
  BASE_::output(os);
  os << ")";
}

template<class T>
Evaluatable<T>* Block<T>::partial_evaluate_replace(const Evaluator<T>& p, bool isarg)
{
  Expression<T>::partial_evaluate(p,isarg);
  return this;
}

} // end namespace expression
} // end namespace alps

#endif // ! ALPS_EXPRESSION_IMPL_H
