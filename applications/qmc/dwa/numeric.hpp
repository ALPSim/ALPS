/*****************************************************************************
*
* ALPS Project Applications: Directed Worm Algorithm for Boson Hubbard Model 
*
* Copyright (C) 2011 by Lode Pollet      <pollet@itp.phys.ethz.ch>    (Inventor), 
*                       Ping Nang Ma     <pingnang@itp.phys.ethz.ch>  (Coder)   ,
*                       Matthias Troyer  <troyer@itp.phys.ethz.ch>    (Advisor) 
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

#ifndef ALPS_APPLICATIONS_NUMERIC_HPP
#define ALPS_APPLICATIONS_NUMERIC_HPP

#include <cassert>
#include <cstdio>
#include <cstdarg>
#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <numeric>
#include <iterator>

#include <boost/bind.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/math/special_functions/round.hpp>


namespace alps {
namespace applications {
namespace numeric {

template <class T>
inline std::vector<T> make_vector(T const x_)  
{  return std::vector<T>(1,x_); }
template <class T>
inline std::vector<T> make_vector(T const x_, T const y_)    
{  std::vector<T> _vec;  _vec.reserve(2);  _vec.push_back(x_);  _vec.push_back(y_);  return _vec;  }
template <class T>
inline std::vector<T> make_vector(T const x_, T const y_, T const z_)   
{  std::vector<T> _vec;  _vec.reserve(3);  _vec.push_back(x_);  _vec.push_back(y_);  _vec.push_back(z_); return _vec;  }

template <class T>
inline T inner_product(std::vector<T> const & vec1_, std::vector<T> const & vec2_)
{  return std::inner_product(vec1_.begin(), vec1_.end(), vec2_.begin(), T());  }

template <class T>
inline T norm(std::vector<T> const & vec_) 
{  return inner_product<T>(vec_,vec_);  }

template <class outType, class inType>
inline std::vector<outType> vector_cast(std::vector<inType> const & x_)  
{  std::vector<outType> y;  y.reserve(x_.size());  std::copy(x_.begin(), x_.end(), std::back_inserter(y));  return y;  }

int pow(int x, int n)  
{ 
  return std::accumulate(boost::counting_iterator<int>(0), boost::counting_iterator<int>(n), 1, boost::bind(std::multiplies<int>(), boost::lambda::_1, x)); 
}
  
inline int iround(double x_)  {  return boost::math::iround(x_);  }

inline double round(double x_, const int places_)
{
  int magnifier = alps::applications::numeric::pow(10,places_);
  int x         = x_ * magnifier;
  return (static_cast<double>(iround(x)) / magnifier);
}

inline std::vector<int>  iround(std::vector<double> x_)  
{  
  std::vector<int> y; 
  y.reserve(x_.size());  
  std::transform(x_.begin(), x_.end(), std::back_inserter(y), static_cast<int (*)(double)>(iround));
  return y;
}

inline std::vector<double> round(std::vector<double> x_, const int places_)
{
  std::vector<double> y; 
  y.reserve(x_.size());
  std::transform
    ( x_.begin(), x_.end(), std::back_inserter(y)
    , boost::bind
        ( static_cast<double (*)(double, const int)>(&round)
        , boost::lambda::_1
        , places_
        )
    );
  return y;
}

template <class T>
inline T mod(T x_)
{
  while(x_<0. || x_>=1.)  
    x_ += (x_>=1. ? -std::floor(x_) : (x_<0. ? -std::ceil(x_)+1 : 0.));
  return x_;
}

} // ending namespace numeric
} // ending namespace applications
} // ending namespace alps

#endif
