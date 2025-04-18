/****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Ping Nang Ma <pingnang@itp.phys.ethz.ch>,
*                            Matthias Troyer <troyer@itp.phys.ethz.ch>,
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

/* $Id: nobinning.h 3520 2009-12-11 16:49:53Z gamperl $ */

#include <alps/numeric/vector_functions.hpp>
#include <iostream>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>


int main(int argc, char** argv)
{
  using namespace alps::numeric;

  // (+,-,*,/) operator(s) for vector are being tested
  std::vector<double> vecA(10,0.81);
  std::vector<double> vecB(10,1.21);
  std::vector<double> res;

  std::cout << "\nA: \n";
  std::copy(vecA.begin(),vecA.end(),std::ostream_iterator<double>(std::cout,"\n"));

  std::cout << "\nB: \n";
  std::copy(vecB.begin(),vecB.end(),std::ostream_iterator<double>(std::cout,"\n"));

  // ( + , - , * , / ) operators for vector are okay
  res = +vecA;
  std::cout << "\n+A: \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = -vecA;
  std::cout << "\n-A: \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = vecA + vecB;
  std::cout << "\nA + B: \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = vecA - vecB;
  std::cout << "\nA - B: \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = vecA * vecB;
  std::cout << "\nA * B: \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = vecA / vecB;
  std::cout << "\nA / B: \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = vecA + 1.;
  std::cout << "\nA + 1.: \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = vecA - 1.;
  std::cout << "\nA - 1.: \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = vecA * 1.;
  std::cout << "\nA * 1.: \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = vecA / 1.;
  std::cout << "\nA / 1.: \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = 1. + vecA;
  std::cout << "\n1. + A: \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = 1. - vecA;
  std::cout << "\n1. - A: \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = 1. * vecA;
  std::cout << "\n1. * A: \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = 1. / vecA;
  std::cout << "\n1. / A: \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));


  // (pow, sq, sqrt, cb, cbrt, exp, log) function(s) for vector are okay
  res = pow(vecA,2.71);
  std::cout << "\nA^(2.71): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = sq(vecA);
  std::cout << "\nsq(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = sqrt(vecA);
  std::cout << "\nsqrt(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = cb(vecA);
  std::cout << "\ncb(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = cbrt(vecA);
  std::cout << "\ncbrt(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = exp(vecA);
  std::cout << "\nexp(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = log(vecA);
  std::cout << "\nlog(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));


  // (sin, cos, tan, asin, acos, atan) function(s) for vector are okay
  res = sin(vecA);
  std::cout << "\nsin(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = cos(vecA);
  std::cout << "\ncos(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = tan(vecA);
  std::cout << "\ntan(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = asin(vecA);
  std::cout << "\nasin(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = acos(vecA);
  std::cout << "\nacos(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = atan(vecA);
  std::cout << "\natan(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));


  // (sinh, cosh, tanh, asinh, acosh, atanh) function(s) for vector are okay
  res = sinh(vecA);
  std::cout << "\nsinh(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = cosh(vecA);
  std::cout << "\ncosh(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = tanh(vecA);
  std::cout << "\ntanh(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = asinh(vecA);
  std::cout << "\nasinh(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = acosh(vecB);
  std::cout << "\nacosh(B): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));

  res = atanh(vecA);
  std::cout << "\natanh(A): \n";
  std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));


  // compound statements
  res = sqrt(vecA*vecA + vecB*vecB);
  std::cout << "\nsqrt(A^2 + B^2): \n";
  //std::copy(res.begin(),res.end(),std::ostream_iterator<double>(std::cout,"\n"));
  std::cout << res << std::endl;

  return 0;
}

