#include "externalsolver.h"

ExternalSolver::multiple_vector_type ExternalSolver::solve(const multiple_vector_type& G0, const alps::Parameters& parms) 
{
  std::string infile = std::tmpnam("impurity_input");
  std::string outfile = std::tmpnam("impurity_output");
 
.... write XML files... into file
  
  system(exe_.native_file_string() + " " + infile + " " + outfile);
 
  multiple_vector_type g;
  .... read g from output file
  
  return g;
}
