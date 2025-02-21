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



#ifndef ALPS_NUMERIC_FOURIER_HPP
#define ALPS_NUMERIC_FOURIER_HPP

#include <alps/numeric/vector_functions.hpp>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>


namespace alps {
namespace numeric {

struct fourier_real
{
  template <class S, class T>
  static T
    evaluate
      (  std::vector<T>  const &  A    //  [a0, a1, a2,...]
      ,  std::vector<T>  const &  B    //  [b0, b1, b2,...]
      ,  T                const &  x    //  x
      ,  S                const & N    //  N
      )
    {
      T  _2PIx_over_N  =  ((M_2PI * x) / N);

      std::vector<T>  _2PInx_over_N;
      _2PInx_over_N.reserve(A.size());
      std::transform
        (  boost::counting_iterator<T>(0)
        ,  boost::counting_iterator<T>(A.size())
        ,  std::back_inserter(_2PInx_over_N)
        , boost::lambda::_1 * _2PIx_over_N
        );  

      std::vector<T>  COS = alps::numeric::cos(_2PInx_over_N);
      std::vector<T>  SIN = alps::numeric::sin(_2PInx_over_N);

      using alps::numeric::operator*;

      COS = COS * A;
      SIN = SIN * B;

      return ( std::accumulate(COS.begin(), COS.end(), 0.) +  std::accumulate(SIN.begin()+1, SIN.end(), 0.) );
    } // objective:  F(x) = a0 + a1 * cos(2\pix/N) + a2 * cos(4\pix/N) + ... + b1 * sin(2\pix/N) + b2 * sin(4\pix/N) + ...
};


} // ending namespace numeric
} // ending namespace alps


#endif

