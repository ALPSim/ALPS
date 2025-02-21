/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2003-2006 by Matthias Troyer <troyer@comp-phys.org>,
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

#ifndef ALPS_MODEL_OPERATOR_H
#define ALPS_MODEL_OPERATOR_H

#include <alps/expression.h>
#include <alps/parameter.h>

namespace alps {

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
