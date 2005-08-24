#include "selfconsistency.h"

void selfconsistency_loop(alps::Parameters& parms, ImpuritySolver& solver, const HilbertTransformer& hilbert)
{
  ImpuritySolver::multiple_vector_type G0,G,Sigma;
  
  G0 = .... initial guess ....;
  
  do {
    g=solver.solve(G0,parms);
    ... hilbert(Sigma[..],mu).....
    ... self consistency step ...
  } while (  ... not converged ....)
}


