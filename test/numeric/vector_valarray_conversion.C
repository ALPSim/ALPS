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

#include <alps/numeric/vector_valarray_conversion.hpp>
#include <iostream>
#include <iterator>
#include <iomanip>
#include <valarray>
#include <vector>
#include <algorithm>


int main(int argc, char** argv)
{
  using namespace alps::numeric;

  std::valarray<double> valA;
  std::vector<double>   vecA;

  valA.resize(10);
  for (int i=0; i < 10; ++i)  {  valA[i] = i;  }

  vecA = valarray2vector<double>(valA);

  std::cout << "Valarray -> Vector\n";
  std::cout << "------------------\n";

  std::cout << "valA:\t";
  std::copy(&valA[0],&valA[valA.size()],std::ostream_iterator<double>(std::cout,"\t"));
  std::cout << "\n";

  std::cout << "vecA:\t";
  std::copy(vecA.begin(),vecA.end(),std::ostream_iterator<double>(std::cout,"\t"));
  std::cout << "\n";


  std::vector<double>   vecB;

  for (int i=0; i < 10; ++i)  {  vecB.push_back(10-i);  }
  
  std::valarray<double> valB = vector2valarray<double>(vecB);

  std::cout << "Vector -> Valarray\n";
  std::cout << "------------------\n";

  std::cout << "vecB:\t";
  std::copy(vecB.begin(),vecB.end(),std::ostream_iterator<double>(std::cout,"\t"));
  std::cout << "\n";

  std::cout << "valB:\t";
  std::copy(&valB[0],&valB[valB.size()],std::ostream_iterator<double>(std::cout,"\t"));
  std::cout << "\n";



  return 0;
}
