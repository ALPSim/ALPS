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
  
  // equilibration
  for (int i=0;i<M/5;++i) {
    if (rng() <= n )
     --n;
    else
      ++n;
  }
  
  histogram.reset(true); 
  
  for (int i=0;i<M;++i) {
    if (rng() <= n )
     --n;
    else
      ++n;

    current_distribution[n]=1;
    histogram << current_distribution;
    current_distribution[n]=0;  
  }
  
  std::cout << "Result: " << histogram << "\n";

  for (int i=0;i<=N;++i)
    std::cout << i << "\t" 
              << probability(N,i) << "\t" 
              << histogram.mean()[i] << "\t" 
              << histogram.error()[i] << "\t" 
              << histogram.tau()[i] << "\n";
  
  return 0;
}
