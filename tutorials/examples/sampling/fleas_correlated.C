/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2007 - 2010 by Matthias Troyer <troyer@comp-phys.org>
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

#include "fleas.h"
#include <alps/alea/observables.hpp>
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
  alps::alea::accumulator histogram(alps::alea::observable("Distribution"));
  alps::alea::accumulator number(alps::alea::observable("Mean number on Anik"));
  
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
  
  const alps::alea::observable_result distribution = histogram.result();
  const std::vector<double> mean = distribution.mean_vector(), error = distribution.error_vector(),
                            tau = distribution.tau_vector();
  std::cout << "Mean number on Anik: " << number.result() << "\n";
  std::cout << "Distribution: " << distribution << "\n";

  for (int i=0;i<=N;++i)
    std::cout << i << "\t" 
              << probability(N,i) << "\t" 
              << mean[i] << "\t" 
              << error[i] << "\t" 
              << tau[i] << "\n";
  
  return 0;
}
