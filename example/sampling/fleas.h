#include <cmath>
double factorial(unsigned int x)
{
  double res=1.;
  while (x)
    res *=x--;
  return res;
}


double probability(int L, int N)
{
  double probability=factorial(L)/factorial(N)/factorial(L-N);
  return probability * std::pow(0.5,L);
}
