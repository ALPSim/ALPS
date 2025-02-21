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

#ifndef ALPS_EXPRESSION_FUNCTION_H
#define ALPS_EXPRESSION_FUNCTION_H

#include <alps/expression/expression_fwd.h>

namespace alps {
namespace expression {


template<class T>
class Function : public Evaluatable<T> {
public:
  typedef T value_type;

  Function(std::istream&, const std::string&);
  Function(const std::string& n, const Expression<T>& e) : name_(n), args_(1,e) {}
  Function(const std::string& n, const std::vector<Expression<T> >& e) : name_(n), args_(e) {}
  value_type value(const Evaluator<T>& =Evaluator<T>(), bool=false) const;
  bool can_evaluate(const Evaluator<T>& =Evaluator<T>(), bool=false) const;
  void output(std::ostream&) const;
  Evaluatable<T>* clone() const { return new Function<T>(*this); }
  boost::shared_ptr<Evaluatable<T> > flatten_one();
  Evaluatable<T>* partial_evaluate_replace(const Evaluator<T>& =Evaluator<T>(), bool=false);
  bool depends_on(const std::string& s) const;
private:
 std::string name_;
 std::vector<Expression<T> > args_;
};

//
// implementation of Function<T>
//

template<class T>
Function<T>::Function(std::istream& in,const std::string& name)
  :  name_(name), args_()
{
  char c;
  in >> c;
  if (c!=')') {
    in.putback(c);
    do {
      args_.push_back(Expression<T>(in));
      in >> c;
    } while (c==',');
    if (c!=')')
      boost::throw_exception(std::runtime_error(std::string("received ") + c + " instead of ) at end of function argument list"));
  }
}

template<class T>
bool Function<T>::depends_on(const std::string& s) const {
  if (name_==s) return true;
  for (typename std::vector<Expression<T> >::const_iterator it=args_.begin();it != args_.end();++it)
    if (it->depends_on(s))
      return true;
  return false;
}

template<class T>
boost::shared_ptr<Evaluatable<T> > Function<T>::flatten_one()
{
  for (typename std::vector<Expression<T> >::iterator it=args_.begin();it != args_.end();++it)
    it->flatten();
  return boost::shared_ptr<Expression<T> >();
}

template<class T>
Evaluatable<T>* Function<T>::partial_evaluate_replace(const Evaluator<T>& p, bool isarg)
{
  p.partial_evaluate_expressions(args_,true);
  return new Block<T>(p.partial_evaluate_function(name_,args_,isarg));
}

template<class T>
typename Function<T>::value_type Function<T>::value(const Evaluator<T>& p, bool isarg) const
{
  value_type val=p.evaluate_function(name_,args_,isarg);
  return val;
}

template<class T>
bool Function<T>::can_evaluate(const Evaluator<T>& p, bool isarg) const
{
  return p.can_evaluate_function(name_,args_,isarg);
}

template<class T>
void Function<T>::output(std::ostream& os) const
{
  os << name_ << "(" << write_vector(args_,", ") << ")";
}

} // end namespace expression
} // end namespace alps

#endif // ! ALPS_EXPRESSION_IMPL_H
