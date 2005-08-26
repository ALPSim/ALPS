#include <alps/parameters.h>
#include "hirschfye.h"

int main(int argc, char** argv)
{

  
  std::string infile = argv[1];
  std::string outfile = argv[2];
  
  // read parameters and G0
  
  multiple_vector_type G0;
  
  HirschFyeSolver solver;
  alps::Parameters parms;
  
  multiple_vector_type g = solver.solve(G0,parms);
  
  //... write g into output file

  return 0;
}
