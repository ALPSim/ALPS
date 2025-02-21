/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2007 - 2010 by Matthias Troyer <troyer@comp-phys.org>
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

#include "fleas.h"
#include <boost/random.hpp>
#include <valarray>

int main()
{
  const int N=50; // total number of fleas
  int M;          // number of samples
  std::cin >> M;
  
  typedef boost::mt19937 engine_type;
  typedef boost::uniform_int<> dist_type;
  
  engine_type engine(3242354);
  dist_type dist(0,1);
  boost::variate_generator<engine_type,dist_type> rng(engine,dist);
  
  std::valarray<unsigned long> histogram(N+1);
  
  for (int i=0;i<M;++i) {
    int n=0;
    for (int j=0;j<N;++j)
      if (rng()==0)
        ++n;
    histogram[n]++;
  }
  
  std::valarray<double> mean(N+1);
  std::valarray<double> error(N+1);
  for (int i=0; i<=N ; ++i) {
    mean[i] = static_cast<double>(histogram[i])/M;
    error[i] = std::sqrt((mean[i]-mean[i]*mean[i])/(M-1));
  }
  
  for (int i=0;i<=N;++i)
    std::cout << i << "\t" 
              << probability(N,i) << "\t"
              << mean[i] << "\t" 
              << error[i] << "\n";
  
  return 0;
}
