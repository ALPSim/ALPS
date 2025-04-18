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
#include <alps/alea.h>
#include <boost/random.hpp>


int main()
{
  const int N=50; // total number of fleas
  
  int M; // number of hops
  std::cout << "How many hops? ";
  std::cin >> M;
  unsigned int seed;
  std::cout << "Random number seed? ";
  std::cin >> seed;
  int n=N; // all fleas on left dog
    
  typedef boost::mt19937 engine_type;
  typedef boost::uniform_int<> dist_type;
  
  typedef boost::mt19937 engine_type;
  engine_type engine(seed);
  dist_type dist(1,N);
  boost::variate_generator<engine_type,dist_type> rng(engine,dist);
  
  std::valarray<double> current_distribution(0.,N+1);
  alps::RealVectorObservable histogram;
  alps::RealObservable number;
  
  // equilibration
  for (int i=0;i<M/5;++i) {
    if (rng() <= n )
     --n;
    else
      ++n;
  }
  
  for (int i=0;i<M;++i) {
    if (rng() <= n )
     --n;
    else
      ++n;

    current_distribution[n]=1;
    histogram << current_distribution;
    number << double(n);
    current_distribution[n]=0;  
  }
  
  std::cout << "Mean number on Anik: " << number << "\n";
  std::cout << "Distribution: " << histogram << "\n";

  for (int i=0;i<=N;++i)
    std::cout << i << "\t" 
              << probability(N,i) << "\t" 
              << histogram.mean()[i] << "\t" 
              << histogram.error()[i] << "\t" 
              << histogram.tau()[i] << "\n";
  
  return 0;
}
