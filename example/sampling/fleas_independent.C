#include "fleas.h"
#include <boost/random.hpp>
#include <valarray>

int main()
{
  const int N=50; // total number of fleas
  
  int M; // number of hops
  std::cout << "How many hops? ";
  std::cin >> M;
  int L; // number of simulations
  std::cout << "How many simulations? ";
  std::cin >> L;
  unsigned int seed;
  std::cout << "Random number seed? ";
  std::cin >> seed;
  
  typedef boost::mt19937 engine_type;
  typedef boost::uniform_int<> dist_type;
  
  engine_type engine;
  engine.seed(seed);
  dist_type dist(1,N);
  boost::variate_generator<engine_type,dist_type> rng(engine,dist);
  
  std::valarray<unsigned long> histogram(N+1);
  std::valarray<double> sum(0.,N+1);
  std::valarray<double> sum_of_squares(0.,N+1);

  // perform L indepedent simulations
  for (int l=0; l<L; ++l) {
  
    int n=N; // all fleas on left dog

    // reset the histogram
    for (int i=0;i <= N; ++i)
      histogram[i]=0.;
      
    // equilibration
    for (int i=0;i<M/5;++i) {
      if (rng() <= n )
       --n;
      else
        ++n;
    }
    
    for (int i=0;i<M/L;++i) {
      if (rng() <= n )
        --n;
      else
        ++n;
      histogram[n]++;  
    }
  
    for (int i=0; i<=N ; ++i) {
      double mean = static_cast<double>(histogram[i])/(M/L);
      sum[i] += mean;
      sum_of_squares[i] += mean*mean;
    }
  }
  
  std::valarray<double> mean(N+1);
  std::valarray<double> error(N+1);
  
  for (int i=0; i<=N ; ++i) {
    mean[i] = sum[i]/L;
    error[i] = std::sqrt((sum_of_squares[i]/L-mean[i]*mean[i])/(L-1));
  }
  
  for (int i=0;i<=N;++i)
    std::cout << i << "\t" 
              << probability(N,i) << "\t"
              << mean[i] << "\t" 
              << error[i] << "\n";
  
  return 0;
}
