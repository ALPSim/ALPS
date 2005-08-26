#ifndef ALPS_DMFT_SOLVER_H
#define ALPS_DMFT_SOLVER_H

#include <alps/parameters.h>
#include <valarray>
#include <utility>

typedef std::valarray<double> vector_type;
typedef std::pair<vector_type,vector_type> multiple_vector_type;

class ImpuritySolver {
public:
  ImpuritySolver() {}
  
  /** solves the impurity problem
  
    @param G0 the inverse Weiss field (without i * omega), as a function of Matsubara frequence omega_n
    @param parms the simulation parameters
    
    Note that we mu, U and h are passed as simulation as simulation parameters.
    
    @return the Green's functions for up and down spins as function of Matsubara frequency */
    
  virtual multiple_vector_type solve(const multiple_vector_type& G0, const alps::Parameters& parms)=0;
  virtual ~ImpuritySolver() {}
};

#endif
