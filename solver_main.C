#include <alps/parameters.h>
#include "selfconsistency.h"

int main(int argc, char** argv)
{

  
  std::string infile = argv[1];
  std::string outfile = argv[2];
  
  read parameters and G0
  
  ImpuritySolver::multiple_vector_type G0;
  
  HirschFyeSolver solver;
  alps::Parameters parms;
  
  ImpuritySolver::multiple_vector_type g = solver.solve(G0,parms);
  
  ... write g into output file

  return 0;
}
