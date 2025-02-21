/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                       Matthias Troyer <troyer@itp.phys.ethz.ch>,
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



#ifndef ALPS_NUMERIC_POLYNOMIAL_HPP
#define ALPS_NUMERIC_POLYNOMIAL_HPP

#include <alps/numeric/vector_functions.hpp>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>


namespace alps {
namespace numeric {


struct polynomial
{
  template <class T>
  static T
    evaluate
      ( std::vector<T>  const &  C
      , T                const &  x
      )
    {
      std::vector<T> X;
      X.reserve(C.size());
      X.push_back(1);
      for (std::size_t index=1; index < C.size(); ++index)
      {
        X.push_back(x * X[index-1]);
      }

      return std::inner_product(C.begin(), C.end(), X.begin(),  0.);
    } // objective: P(x) = C0 + C1 * x + C2 * x^2 + ... 


  template <class T>
  static T
    evaluate_derivative
      ( std::vector<T>  const &  C
      , T                const &  x
      )
    {
      std::vector<T>   X;
      X.reserve(C.size());
      X.push_back(1);
      for (std::size_t index=1; index < C.size()-1; ++index)
      {
        X.push_back(x * X[index-1]);
      }
      
      std::vector<T>  C_prime;
      C_prime.reserve(C.size()-1);
      std::transform(boost::counting_iterator<unsigned int>(1), boost::counting_iterator<unsigned int>(C.size()), C.begin()+1, std::back_inserter(C_prime), boost::lambda::_1 * boost::lambda::_2);

      return std::inner_product(C_prime.begin(), C_prime.end(), X.begin(), 0.);
    }  
};


} // ending namespace numeric
} // ending namespace alps


#endif

