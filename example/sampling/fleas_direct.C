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
