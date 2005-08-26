#include "externalsolver.h"

multiple_vector_type ExternalSolver::solve(const multiple_vector_type& G0, const alps::Parameters& parms) 
{
  std::string infile = "infile"; // = std::tmpnam("infile");
  std::string outfile = "outfile"; // = std::tmpnam("outfile);
 
//.... write XML files... into file
  
  std::string command = exe_.native_file_string();
  std::cerr << command << "\n";
  command += " " + infile + " " + outfile;
  
  system(command.c_str());
 
  multiple_vector_type g;
 // .... read g from output file
  
  return g;
}
